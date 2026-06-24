# This Python file uses the following encoding: utf-8

import psutil
import os
import threading
import multiprocessing
import traceback

if __name__ == '__main__':
    multiprocessing.freeze_support()

################# CLEANUP ####################
"""
Children discovery logic is extracted into a helper function that returns:

Immediate children

Non-immediate (recursive) descendants excluding immediate children

Before each timeout, the list of non-immediate children is re-evaluated.

On the final timeout, force kills all non-immediate children.

If anything is still alive after force kill, it then force kills immediate children as a fallback.
"""

def get_children_immediate_nonimmediate(root_pid = None):
    """
    Returns a tuple of:
    - immediate_children: List[psutil.Process]
    - non_immediate_descendants: List[psutil.Process]
    """
    if root_pid is None:
        root_pid = os.getpid()
    parent_proc = psutil.Process(root_pid)

    immediate = parent_proc.children(recursive=False)
    all_descendants = parent_proc.children(recursive=True)
    immediate_pids = {p.pid for p in immediate}

    non_immediate = [p for p in all_descendants if p.pid not in immediate_pids]
    return immediate, non_immediate


def cleanup_non_immediate_descendants(
      main_pid=None, timeout_sequence=[3,0.25,0.25,0.25,0.25,3]):
    """
    Cleans up non-immediate descendants in stages with increasing timeouts.
    Default timeout sequence: SCIP (still?) needs 5x Ctrl-C.
    Final step: force-kills immediate children if anything is still alive.
    """
    if not timeout_sequence:
        print("[CLEANUP] No timeouts provided. Nothing to do.")
        return

    if main_pid is None:
        main_pid = os.getpid()

    non_immed_0 = set()         ## Linux: need to remember all

    for i, timeout in enumerate(timeout_sequence):
        _, non_immediate = get_children_immediate_nonimmediate(main_pid)
        non_immed_now = non_immediate
        non_immediate = set(non_immediate) | non_immed_0
        if not non_immediate:
            print(f"[CLEANUP] No non-immediate child processes at stage {i+1}.")
            return

        non_immed_0 = non_immediate
        is_final = i == len(timeout_sequence) - 1

        action = "KILL" if is_final else "TERMINATE"
        print(f"[CLEANUP] Stage {i+1}: {action} {len(non_immediate)} non-immediate processes (timeout={timeout}s)...")

        killedSome = False
        for proc in non_immediate:
            try:
                print(f"  [{action}] PID={proc.pid}, name={proc.name()}")
                proc.kill() if is_final else proc.terminate()
                killedSome = True
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
            except Exception as e:
                print(f"  [ERROR] Failed to {action.lower()} PID {proc.pid}: {e}")

        if not killedSome and not non_immed_now:
            print(f"[CLEANUP] Stage {i+1} complete. All targeted processes terminated.")
            return

        gone, alive = psutil.wait_procs(non_immediate, timeout=timeout)

        if not alive:
            print(f"[CLEANUP] Stage {i+1} complete. All targeted processes terminated.")
        else:
            print(f"[CLEANUP] Stage {i+1} complete. {len(alive)} processes still alive.")

    # Final check — if anything is still alive, also kill immediate children
    _, remaining_non_immediate = get_children_immediate_nonimmediate(main_pid)
    if remaining_non_immediate:
        print(f"[CLEANUP] Final sweep: {len(remaining_non_immediate)} non-immediate processes still alive.")

        immediate, _ = get_children_immediate_nonimmediate(main_pid)
        if immediate:
            print(f"[CLEANUP] Killing {len(immediate)} immediate child processes as fallback.")
            for proc in immediate:
                try:
                    print(f"  [KILL IMMEDIATE] PID={proc.pid}, name={proc.name()}")
                    proc.kill()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    continue
                except Exception as e:
                    print(f"  [ERROR] Failed to kill PID {proc.pid}: {e}")
            psutil.wait_procs(immediate, timeout=3)

    # Final summary
    _, final_non_immediate = get_children_immediate_nonimmediate(main_pid)
    if final_non_immediate:
        print(f"[CLEANUP] WARNING: {len(final_non_immediate)} non-immediate process(es) STILL alive.")
        for proc in final_non_immediate:
            print(f"  [STILL ALIVE] PID={proc.pid}, name={proc.name()}")



###################### Function timeout with watcher in a separate process #########
"""
This is the only design that works reliably in Python 3.13.

threading's timers fail sometimes. Might improve in Python 3.14.

Use multiprocessing.Event for coordination.

Start a watcher process that:

Waits for first_timeout.

If not completed, runs the cleanup.

@todo Can the watcher be a child? But then an immediate one.

Then waits second_timeout, and warns if the target is still running.

"""

def _watcher_process(done_event, cleanup_func, main_pid, first_timeout, second_timeout):
    # Wait up to the first timeout
    if not done_event.wait(timeout=first_timeout):
        print(f"[WATCHER] First timeout ({first_timeout}s) reached. Running cleanup...")
        try:
            cleanup_func(main_pid)
        except Exception as e:
            print(f"[WATCHER] Cleanup function raised exception: {e}")

        # Wait for the second timeout
        if not done_event.wait(timeout=second_timeout):
            print(f"[WATCHER] Second timeout ({second_timeout}s) reached. Target still running.")
        else:
            print("[WATCHER] Target finished during second timeout.")
    else:
        pass
        ## print("[WATCHER] Target finished before first timeout.")

def run_with_timeout_cleanup_multiprocessing_watcher(
      target_func, cleanup_func, first_timeout, second_timeout, *args, **kwargs):
    """
    Runs `target_func` in the main process/thread.
    Starts a multiprocessing watcher that:
      - Waits up to `first_timeout`
      - Runs `cleanup_func` if target not done
      - Waits `second_timeout`
      - Logs if still running
    """
    ctx = multiprocessing.get_context("spawn")  # Safer across platforms
    done_event = ctx.Event()

    main_pid = multiprocessing.current_process().pid
    watcher = ctx.Process(
        target=_watcher_process,
        args=(done_event, cleanup_func, main_pid, first_timeout, second_timeout),
    )
    watcher.start()

    result = None
    try:
        result = target_func(*args, **kwargs)
    except Exception:
        print("[ERROR] Exception occurred in target function:")
        traceback.print_exc()
    finally:
        done_event.set()
        watcher.join()

    return result



###################### Function timeout with watcher in a thread ##############
"""
The timer still unreliable on TotalVariation2D...
"""

def run_with_timeout_cleanup_thread(
      target_func, cleanup_func, first_timeout, second_timeout, *args, **kwargs):
    """
    Runs `target_func` in the main thread.
    Starts a watcher thread that:
      - Waits for `first_timeout`
      - If the target hasn't finished, runs `cleanup_func`
      - Waits for `second_timeout`
      - If the target still hasn't finished, logs it

    If the target finishes early, the watcher exits immediately.
    """
    result_container = {'done': False, 'result': None, 'error': None}
    done_event = threading.Event()

    def watcher():
        # Wait up to first timeout
        if not done_event.wait(timeout=first_timeout):
            print(f"[WATCHER] First timeout ({first_timeout}s) reached. Running cleanup...")
            try:
                cleanup_func()
            except Exception as e:
                print(f"[WATCHER] Cleanup function raised exception: {e}")

            # Wait up to second timeout
            if not done_event.wait(timeout=second_timeout):
                print(f"[WATCHER] Second timeout ({second_timeout}s) reached. Target still running.")
            else:
                print("[WATCHER] Target finished during second timeout.")
        else:
            # Target completed within first timeout
            print("[WATCHER] Target completed before first timeout. Watcher exiting.")

    watch_thread = threading.Thread(target=watcher, daemon=True)
    watch_thread.start()

    try:
        result_container['result'] = target_func(*args, **kwargs)
    except Exception:
        result_container['error'] = traceback.format_exc()
    finally:
        result_container['done'] = True
        done_event.set()

    watch_thread.join()

    if result_container['error']:
        print("[ERROR] Target function raised an exception:")
        print(result_container['error'])

    return result_container['result']



###################### Function timeout with threading ########################
"""
Here’s a threading-based version of the utility. Note that while this works
for monitoring and cleanup, Python threads cannot be forcefully killed, so
if the target function hangs indefinitely (e.g., due to a blocking C extension
or deadlock), the thread will remain alive. Use this version only if you don’t
need to force-terminate the function.

Limitations of the Threading Version:
No way to forcefully kill the thread.

If the thread blocks indefinitely, it cannot be recovered or stopped from outside.

Use for non-blocking, cooperative code only (e.g., functions you expect to be
well-behaved or time-limited).

If you ever need robust cancellation, stick with the multiprocessing version (below).
"""

def run_with_timeout_threaded(target_func, cleanup_func, first_timeout, second_timeout, *args, **kwargs):
    """
    Runs `target_func` in a thread. If it doesn't complete within `first_timeout` seconds,
    runs `cleanup_func`, then waits another `second_timeout` seconds.
    Cannot force-terminate threads, so this will only monitor and report.
    """
    result_container = {'result': None, 'error': None}

    def thread_wrapper():
        try:
            result_container['result'] = target_func(*args, **kwargs)
        except Exception:
            result_container['error'] = traceback.format_exc()

    thread = threading.Thread(target=thread_wrapper)
    print("START TARGET WITH tm", first_timeout)
    thread.start()
    print("TARGET STARTED. JOINING")

    thread.join(timeout=first_timeout)
    print("JOINED TARGET WITH tm", first_timeout)

    if thread.is_alive():
        print(f"[INFO] Function did not complete in {first_timeout}s. Running cleanup...")
        try:
            cleanup_func()
        except Exception as e:
            print(f"[WARNING] Cleanup function raised an exception: {e}")

        thread.join(timeout=second_timeout)

        if thread.is_alive():
            print(f"[WARNING] Function still running after additional {second_timeout}s. Cannot terminate thread.")
            return None

    if result_container['error']:
        print("[ERROR] Function raised an exception:")
        print(result_container['error'])
        return None

    return result_container['result']


###################### Function timeout with multiprocessing ##################
"""
Here's a cross-platform Python utility that:

Executes a given function in a separate thread or process.

Waits for a first timeout.

If the function hasn’t completed by then, it executes a cleanup function.

Then waits again for a second timeout for the original function to finish.

If the function still doesn't finish, it terminates the process.

This approach uses the multiprocessing module for cross-platform compatibility
(works on Windows and Unix).

Notes:
This utility ensures true isolation of the target function using multiprocessing
so it can be terminated cleanly on any OS.

Threads don’t support reliable termination in Python, so multiprocessing is used instead.

Cleanup happens in the main process, not inside the child, ensuring access to resources.

"""

def wrapper(target_func, queue, *args, **kwargs):
    try:
        result = target_func(*args, **kwargs)
        queue.put(('result', result))
    except Exception as e:
        queue.put(('error', traceback.format_exc()))


def run_with_timeout_multiprocessing(
      target_func, cleanup_func, first_timeout, second_timeout, *args, **kwargs):
    """
    Runs `target_func` with given arguments. If it doesn't complete within `first_timeout` seconds,
    calls `cleanup_func`, then waits another `second_timeout` seconds before force-terminating.
    """
    queue = multiprocessing.Queue()
    process = multiprocessing.Process(target=wrapper,
        args=(target_func, queue, *args), kwargs=kwargs)
    process.start()

    process.join(timeout=first_timeout)

    if process.is_alive():
        print(f"[INFO] Function did not complete in {first_timeout}s. Running cleanup...")
        try:
            cleanup_func()
        except Exception as e:
            print(f"[WARNING] Cleanup function raised an exception: {e}")

        process.join(timeout=second_timeout)

        if process.is_alive():
            print(f"[WARNING] Function still did not finish after additional {second_timeout}s. Terminating...")
            process.terminate()
            process.join()
            return None

    if not queue.empty():
        status, data = queue.get()
        if status == 'result':
            return data
        elif status == 'error':
            print("[ERROR] Function raised an exception:")
            print(data)
            return None

    print("[INFO] Function finished but returned no result.")
    return None
