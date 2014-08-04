/*
 Classes for interacting with Java code through JNI.

 Copyright (C) 2012 AMPL Optimization Inc

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization Inc disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "solvers/jacop/java.h"

#include <cstdlib>
#include <vector>

#ifdef _WIN32
# include <windows.h>
# include <sys/stat.h>
#endif

namespace {

class String {
 private:
  JNIEnv *env_;
  jstring str_;
  const char *utf_chars_;

  FMT_DISALLOW_COPY_AND_ASSIGN(String);

 public:
  String(JNIEnv *env, jstring s) :
    env_(env), str_(s), utf_chars_(env->GetStringUTFChars(s, 0)) {
  }

  ~String() {
    env_->ReleaseStringUTFChars(str_, utf_chars_);
  }

  const char *c_str() const {
    return utf_chars_;
  }
};

#ifdef _WIN32
// A registry key.
class RegKey {
 private:
  HKEY key_;
  
  FMT_DISALLOW_COPY_AND_ASSIGN(RegKey);

 public:
  RegKey(HKEY key, fmt::StringRef subkey, REGSAM access);
  ~RegKey();

  HKEY get() const { return key_; }

  std::string GetSubKeyName(int index) const;
  std::string GetStrValue(fmt::StringRef name) const;
};

RegKey::RegKey(HKEY key, fmt::StringRef subkey, REGSAM access) : key_() {
  LONG result = RegOpenKeyExA(key, subkey.c_str(), 0, access, &key_);
  if (result != ERROR_SUCCESS) {
    throw fmt::WindowsError(GetLastError(),
      "cannot open registry key {}", subkey);
  }
}

RegKey::~RegKey() { RegCloseKey(key_); }

std::string RegKey::GetSubKeyName(int index) const {
  DWORD max_subkey_len = 0;
  LONG result = RegQueryInfoKeyA(key_,
      0, 0, 0, 0, &max_subkey_len, 0, 0, 0, 0, 0, 0);
  if (result != ERROR_SUCCESS)
    max_subkey_len = 256;
  std::vector<char> name(max_subkey_len + 1);
  DWORD name_size = static_cast<DWORD>(name.size());
  result = RegEnumKeyExA(key_, 0, &name[0], &name_size, 0, 0, 0, 0);
  if (result != ERROR_SUCCESS)
    throw fmt::WindowsError(GetLastError(), "cannot get registry key name");
  return &name[0];
}

std::string RegKey::GetStrValue(fmt::StringRef name) const {
  char buffer[256];
  DWORD size = sizeof(buffer);
  DWORD type = 0;
  LONG result = RegQueryValueExA(key_,
      name.c_str(), 0, &type, reinterpret_cast<LPBYTE>(buffer), &size);
  if (result != ERROR_SUCCESS) {
    throw fmt::WindowsError(GetLastError(),
      "cannot get registry key value {}", name);
  }
  if (type != REG_SZ) {
    throw fmt::WindowsError(
      GetLastError(), "value of key {} is not a string", name);
  }
  buffer[std::min<DWORD>(size, sizeof(buffer) - 1)] = 0;
  return buffer;
}
#endif
}

namespace ampl {

JVM JVM::instance_;

void Env::Throw(jthrowable exception, const char *method_name) {
  env_->ExceptionClear();
  jmethodID getMessage = GetMethod(FindClass("java/lang/Object"),
      "toString", "()Ljava/lang/String;");
  String message(env_, static_cast<jstring>(Check(
      env_->CallObjectMethod(exception, getMessage), "CallObjectMethod")));
  throw JavaError(fmt::format(
      "{} failed: {}", method_name, message.c_str()), exception);
}

jobject Env::NewObject(jclass cls, jmethodID ctor, ...) {
  std::va_list args;
  va_start(args, ctor);
  jobject result = NewObjectV(cls, ctor, args);
  va_end(args);
  return Check(result, "NewObjectV");
}

jobject Env::NewObject(const char *class_name, const char *ctor_sig, ...) {
  jclass cls = FindClass(class_name);
  jmethodID ctor = GetMethod(cls, "<init>", ctor_sig);
  std::va_list args;
  va_start(args, ctor_sig);
  jobject result = env_->NewObjectV(cls, ctor, args);
  va_end(args);
  return Check(result, "NewObjectV");
}

void Env::CallVoidMethod(jobject obj, jmethodID method, ...) {
  std::va_list args;
  va_start(args, method);
  env_->CallVoidMethodV(obj, method, args);
  va_end(args);
  Check("CallVoidMethodV");
}

jboolean Env::CallBooleanMethod(jobject obj, jmethodID method, ...) {
  std::va_list args;
  va_start(args, method);
  jboolean result = env_->CallBooleanMethodV(obj, method, args);
  va_end(args);
  Check("CallBooleanMethodV");
  return result;
}

jint Env::CallIntMethod(jobject obj, jmethodID method, ...) {
  std::va_list args;
  va_start(args, method);
  jint result = env_->CallIntMethodV(obj, method, args);
  va_end(args);
  Check("CallIntMethodV");
  return result;
}

jint Env::CallIntMethodKeepException(jobject obj, jmethodID method, ...) {
  std::va_list args;
  va_start(args, method);
  jint result = env_->CallIntMethodV(obj, method, args);
  va_end(args);
  if (env_->ExceptionOccurred())
    throw JavaError("CallIntMethodV failed");
  return result;
}

JVM::~JVM() {
  if (jvm_)
    jvm_->DestroyJavaVM();
}

Env JVM::env(const char *const *options) {
  if (!instance_.jvm_) {
#ifdef _WIN32
    std::string runtime_lib_path;
    bool exists = false;
    try {
      RegKey jre_key(HKEY_LOCAL_MACHINE,
          "SOFTWARE\\JavaSoft\\Java Runtime Environment",
          KEY_ENUMERATE_SUB_KEYS);
        RegKey key(jre_key.get(), jre_key.GetSubKeyName(0), KEY_QUERY_VALUE);
        runtime_lib_path = key.GetStrValue("RuntimeLib");
        struct _stat s = {};
        exists = _stat(runtime_lib_path.c_str(), &s) == 0;
    } catch (const JavaError &) {
      // Ignore error.
    }
    std::string::size_type pos = runtime_lib_path.rfind('\\');
    if (pos != std::string::npos) {
      runtime_lib_path = runtime_lib_path.substr(0, pos);
      if (!exists) {
        // A workaround for broken path to jvm.dll on 64-bit Windows.
        pos = runtime_lib_path.rfind('\\');
        if (pos != std::string::npos)
          runtime_lib_path.replace(pos + 1, std::string::npos, "server");
      }
      std::string path = std::getenv("PATH");
      path += ";";
      path += runtime_lib_path;
      path += ";";
      path += runtime_lib_path + "\\..";
      SetEnvironmentVariable("PATH", path.c_str());
    }
    if (!LoadLibrary("jvm.dll"))
       throw JavaError("Failed to load jvm.dll");
#endif
    JavaVMInitArgs vm_args = {};
    vm_args.version = JNI_VERSION_1_6;
    vm_args.ignoreUnrecognized = false;
    std::vector<JavaVMOption> jvm_options;
    if (options) {
      for (const char *const *opt = options; *opt; ++opt) {
        JavaVMOption jvm_opt = {const_cast<char*>(*opt)};
        jvm_options.push_back(jvm_opt);
      }
    }
    vm_args.nOptions = static_cast<jint>(jvm_options.size());
    vm_args.options = &jvm_options[0];
    void *envp = 0;
    jint result = JNI_CreateJavaVM(&instance_.jvm_, &envp, &vm_args);
    if (result != JNI_OK) {
      throw JavaError(fmt::format(
          "Java VM initialization failed, error code = {}", result));
    }
    instance_.env_ = Env(static_cast<JNIEnv*>(envp));
  }
  return instance_.env_;
}

ClassBase::~ClassBase() {}

jobject ClassBase::NewObject(Env env, ...) {
  Init(env);
  std::va_list args;
  va_start(args, env);
  jobject result = 0;
  try {
    result = env.NewObjectV(class_, ctor_, args);
    va_end(args);
  } catch (...) {
    va_end(args);
    throw;
  }
  return result;
}
}
