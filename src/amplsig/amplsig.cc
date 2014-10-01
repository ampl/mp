/*
 AMPL extension DLL that simulates sending SIGINT (Ctrl-C) from another
 process on Windows. It is required since GenerateConsoleCtrlEvent doesn't
 work reliably with CTRL_C_EVENT across processes.
 The DLL can be loaded from AMPL as follows:
   load signal.dll;
 It connects to a socket at 127.0.0.1:$AMPLSIG_PORT and generates
 CTRL_C_EVENT on every SIGINT command read from the socket input stream.

 Copyright (C) 2014 AMPL Optimization Inc

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

#include <cstdio>
#include <cstdlib>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>

#pragma comment(lib, "Ws2_32.lib")

using namespace std;

BOOL APIENTRY DllMain(HMODULE, DWORD, LPVOID) {
  return TRUE;
}

class Thread {
 private:
  WSADATA wsa_data_;
  addrinfo* addr_;
  SOCKET socket_;
  int debug_level_;

  // Do not implement!
  Thread(const Thread&);
  void operator=(const Thread&);

  void Fail(const char* func_name, int code);
  void Log(const char* message) {
    if (debug_level_ > 1)
      fprintf(stderr, "amplsig: %s\n", message);
  }

  Thread();
  ~Thread();

  void DoRun();

 public:
  static DWORD WINAPI Run(LPVOID);
};

void Thread::Fail(const char* func_name, int code) {
  if (debug_level_ > 0)
    fprintf(stderr, "amplsig: %s failed: %d\n", func_name, code);
  throw code;
}

Thread::Thread() : addr_(0), socket_(INVALID_SOCKET), debug_level_(0) {
  const char* debug = getenv("AMPLSIG_DEBUG");
  if (debug) debug_level_ = atoi(debug);
  int result = WSAStartup(MAKEWORD(2, 2), &wsa_data_);
  if (result != 0)
    Fail("WSAStartup", result);
}

Thread::~Thread() {
  if (socket_ != INVALID_SOCKET)
    closesocket(socket_);
  if (addr_)
    freeaddrinfo(addr_);
  WSACleanup();
}

void Thread::DoRun() {
  const char* port = getenv("AMPLSIG_PORT");
  if (!port) {
    Log("Port is not specified");
    return;
  }

  addrinfo hints = {0};
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  hints.ai_protocol = IPPROTO_TCP;
  int result = getaddrinfo("127.0.0.1", port, &hints, &addr_);
  if (result != 0)
    Fail("getaddrinfo", result);
  socket_ = ::socket(addr_->ai_family,
    addr_->ai_socktype, addr_->ai_protocol);
  if (socket_ == INVALID_SOCKET)
    Fail("socket", WSAGetLastError());
  Log("Connecting");
  result = connect(socket_, addr_->ai_addr, (int)addr_->ai_addrlen);
  if (result == SOCKET_ERROR)
    Fail("connect", WSAGetLastError());
  Log("Connected");

  enum { BUFFER_LEN = 512 };
  char buffer[BUFFER_LEN];
  do {
    result = recv(socket_, buffer, BUFFER_LEN, 0);
    if (result > 0) {
      if (strncmp(buffer, "SIGINT", 6) == 0 &&
          !GenerateConsoleCtrlEvent(CTRL_C_EVENT, 0)) {
        Fail("GenerateConsoleCtrlEvent", GetLastError());
      }
    } else if (result != 0)
      Fail("recv", WSAGetLastError());
  } while (result > 0);
}

DWORD WINAPI Thread::Run(LPVOID) {
  try {
    Thread().DoRun();
  } catch (...) {}
  return 0;
}

extern "C" __declspec(dllexport) void funcadd_ASL(void*) {
  CreateThread(NULL, 0, Thread::Run, NULL, 0, NULL);
}
