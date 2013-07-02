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

#ifndef AMPL_SOLVERS_JACOP_JAVA_H_
#define AMPL_SOLVERS_JACOP_JAVA_H_

#include <jni.h>
#include <cstdarg>
#include <stdexcept>
#include <string>

#include "solvers/util/format.h"

#ifdef WIN32
# define AMPL_CLASSPATH_SEP ";"
#else
# define AMPL_CLASSPATH_SEP ":"
#endif

namespace ampl {

class JavaError : public std::runtime_error {
 private:
  jthrowable exception_;

 public:
  explicit JavaError(fmt::StringRef message, jthrowable exception = 0)
    : std::runtime_error(message), exception_(exception) {}

  jthrowable exception() const { return exception_; }
};

// Java Native Interface environment.
class Env {
 private:
  JNIEnv *env_;

  void Throw(jthrowable exception, const char *method_name);

  // Checks the result of a method call returning pointer and throws an
  // exception in the case of a failure.
  template<typename T>
  T Check(T result, const char *method_name);

  void Check(const char *method_name) {
    if (jthrowable exception = env_->ExceptionOccurred())
      Throw(exception, method_name);
  }

 public:
  explicit Env(JNIEnv *env = 0) : env_(env) {}

  jclass FindClass(const char *name) {
    return Check(env_->FindClass(name), "FindClass");
  }

  jmethodID GetMethod(jclass cls, const char *name, const char *sig) {
    return Check(env_->GetMethodID(cls, name, sig), "GetMethodID");
  }

  jfieldID GetFieldID(jclass cls, const char *name, const char *sig) {
    return Check(env_->GetFieldID(cls, name, sig), "GetFieldID");
  }

  jfieldID GetStaticFieldID(jclass cls, const char *name, const char *sig) {
    return Check(env_->GetStaticFieldID(cls, name, sig), "GetStaticFieldID");
  }

  jboolean GetBooleanField(jobject obj, jfieldID field) {
    jboolean value = env_->GetBooleanField(obj, field);
    Check("GetBooleanField");
    return value;
  }

  jint GetStaticIntField(jclass cls, jfieldID field) {
    jint value = env_->GetStaticIntField(cls, field);
    Check("GetStaticIntField");
    return value;
  }

  jobject NewObject(jclass cls, jmethodID ctor, ...);

  jobject NewObjectV(jclass cls, jmethodID ctor, std::va_list args) {
    return Check(env_->NewObjectV(cls, ctor, args), "NewObjectV");
  }

  jobject NewObject(const char *class_name, const char *ctor_sig, ...);

  void CallVoidMethod(jobject obj, jmethodID method, ...);
  jboolean CallBooleanMethod(jobject obj, jmethodID method, ...);
  jint CallIntMethod(jobject obj, jmethodID method, ...);

  jobjectArray NewObjectArray(jsize length,
      jclass elementClass, jobject initialElement) {
    return Check(env_->NewObjectArray(length, elementClass, initialElement),
        "NewObjectArray");
  }

  void SetObjectArrayElement(jobjectArray array, jsize index, jobject value) {
    env_->SetObjectArrayElement(array, index, value);
    Check("SetObjectArrayElement");
  }

  jintArray NewIntArray(jsize length) {
    return Check(env_->NewIntArray(length), "NewIntArray");
  }

  void SetIntArrayRegion(jintArray array,
      jsize start, jsize length, const jint *values) {
    env_->SetIntArrayRegion(array, start, length, values);
    Check("SetIntArrayRegion");
  }

  jboolean IsInstanceOf(jobject obj, jclass cls) {
    jboolean result = env_->IsInstanceOf(obj, cls);
    Check("IsInstanceOf");
    return result;
  }

  void RegisterNatives(jclass cls, const JNINativeMethod *methods, jint size) {
    if (env_->RegisterNatives(cls, methods, size) < 0)
      Check(0, "RegisterNatives");
  }
};

template <typename T>
T Env::Check(T result, const char *method_name) {
  if (!result) {
    jthrowable exception = env_->ExceptionOccurred();
    if (!exception)
      throw JavaError(std::string(method_name) + " failed");
    Throw(exception, method_name);
  }
  return result;
}

// Java Virtual Machine.
// A process can have at most one JVM, therefore there is a single static JVM
// object which is initialized when JVM::env() is called for the first time.
class JVM {
 private:
  JavaVM *jvm_;
  Env env_;
  static JVM instance_;

  JVM() : jvm_() {}
  ~JVM();

 public:
  // options: an array of JVM options terminated by a null pointer.
  static Env env(const char *const *options = 0);
};

class ClassBase {
 protected:
  jclass class_;
  jmethodID ctor_;

  void True() const {}
  typedef void (ClassBase::*SafeBool)() const;

  virtual void DoInit(Env env) = 0;

 public:
  ClassBase() : class_(), ctor_() {}
  virtual ~ClassBase();

  void Init(Env env) {
    if (!class_)
      DoInit(env);
  }

  operator SafeBool() const { return class_ ? &ClassBase::True : 0; }

  jclass get() const { return class_; }

  jobject NewObject(Env env, ...);
};

// A reference to a Java class and one of its constructor.
template <typename Info>
class Class : public ClassBase {
 protected:
  void DoInit(Env env) {
    class_ = env.FindClass(Info::name());
    ctor_ = env.GetMethod(class_, "<init>", Info::ctor_sig());
  }
};

#define CLASS_INFO(class_name, jvm_class_name, jvm_ctor_sig) \
struct class_name { \
  static const char *name() { return jvm_class_name; } \
  static const char *ctor_sig() { return jvm_ctor_sig; } \
};
}

#endif // AMPL_SOLVERS_JACOP_JAVA_H_
