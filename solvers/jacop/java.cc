/*
 Classes for interacting with Java code through JNI.

 Copyright (C) 2012 AMPL Optimization LLC

 Permission to use, copy, modify, and distribute this software and its
 documentation for any purpose and without fee is hereby granted,
 provided that the above copyright notice appear in all copies and that
 both that the copyright notice and this permission notice and warranty
 disclaimer appear in supporting documentation.

 The author and AMPL Optimization LLC disclaim all warranties with
 regard to this software, including all implied warranties of
 merchantability and fitness.  In no event shall the author be liable
 for any special, indirect or consequential damages or any damages
 whatsoever resulting from loss of use, data or profits, whether in an
 action of contract, negligence or other tortious action, arising out
 of or in connection with the use or performance of this software.

 Author: Victor Zverovich
 */

#include "solvers/jacop/java.h"

#include "solvers/util/noncopyable.h"

#include <cstdlib>
#include <vector>

#ifdef WIN32
# include <windows.h>
# define CLASSPATH_SEP ";"
#else
# define CLASSPATH_SEP ":"
#endif

namespace {

class String : ampl::Noncopyable {
 private:
  JNIEnv *env_;
  jstring str_;
  const char *utf_chars_;

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

#ifdef WIN32
// A registry key.
class RegKey : ampl::Noncopyable {
 private:
  HKEY key_;

 public:
  RegKey(HKEY key, const char *subkey, REGSAM access);
  ~RegKey();

  std::string GetSubKeyName(int index) const;
  std::string GetStrValue(fmt::StringRef subkey, fmt::StringRef name) const;
};

RegKey::RegKey(HKEY key, const char *subkey, REGSAM access) : key_() {
  LONG result = RegOpenKeyExA(key, subkey, 0, access, &key_);
  if (result != ERROR_SUCCESS) {
    throw ampl::JavaError(fmt::Format(
        "RegOpenKeyExA failed: error code = {}") << result);
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
  DWORD name_size = name.size();
  result = RegEnumKeyExA(key_, 0, &name[0], &name_size, 0, 0, 0, 0);
  if (result != ERROR_SUCCESS) {
    throw ampl::JavaError(fmt::Format(
        "RegEnumKeyExA failed: error code = {}") << result);
  }
  return &name[0];
}

std::string RegKey::GetStrValue(fmt::StringRef subkey, fmt::StringRef name) const {
  char buffer[256];
  DWORD size = sizeof(buffer);
  LONG result = RegGetValueA(key_,
      subkey.c_str(), name.c_str(), RRF_RT_REG_SZ, 0, buffer, &size);
  if (result != ERROR_SUCCESS) {
    throw ampl::JavaError(fmt::Format(
        "RegGetValueA failed: error code = {}") << result);
  }
  buffer[sizeof(buffer) - 1] = 0;
  return buffer;
}
#endif
}

namespace ampl {

JVM JVM::instance_;

void Env::Throw(jthrowable exception, const char *method_name) {
  jmethodID getMessage = GetMethod(FindClass("java/lang/Throwable"),
      "getMessage", "()Ljava/lang/String;");
  String message(env_, static_cast<jstring>(Check(
      env_->CallObjectMethod(exception, getMessage), "CallObjectMethod")));
  throw JavaError(fmt::Format("{} failed: {}")
        << method_name << message.c_str());
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

JVM::JVM() : jvm_() {
#ifdef WIN32
  std::string runtime_lib_path;
  try {
    RegKey key(HKEY_LOCAL_MACHINE,
        "SOFTWARE\\JavaSoft\\Java Runtime Environment", KEY_ENUMERATE_SUB_KEYS);
    runtime_lib_path = key.GetStrValue(key.GetSubKeyName(0), "RuntimeLib");
  } catch (const JavaError &e) {
    return;  // Ignore error.
  }
  std::string::size_type pos = runtime_lib_path.rfind('\\');
  if (pos == std::string::npos)
    return;
  runtime_lib_path = runtime_lib_path.substr(0, pos);
  std::string path = std::getenv("PATH");
  path += ";";
  path += runtime_lib_path;
  SetEnvironmentVariable("PATH", path.c_str());
#endif
}

JVM::~JVM() {
  if (jvm_)
    jvm_->DestroyJavaVM();
}

Env JVM::env() {
  if (!instance_.jvm_) {
    JavaVMInitArgs vm_args = {};
    vm_args.version = JNI_VERSION_1_6;
    vm_args.ignoreUnrecognized = false;
    JavaVMOption option = {};
    option.optionString = const_cast<char*>(
        "-Djava.class.path=JaCoP-3.2.jar" CLASSPATH_SEP "lib/JaCoP-3.2.jar");
    vm_args.nOptions = 1;
    vm_args.options = &option;
    void *envp = 0;
    jint result = JNI_CreateJavaVM(&instance_.jvm_, &envp, &vm_args);
    if (result != JNI_OK) {
      throw JavaError(fmt::Format(
          "Java VM initialization failed, error code = {}") << result);
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
