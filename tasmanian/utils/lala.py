#!/usr/bin/python

from multiprocessing import Process, Queue
import time, random

def do_something(n_order, x, queue):
    time.sleep(5)
    queue.put((n_order, x))


def main():

    data = [1,2,3,4,5]
    queue = Queue()
    processes = [Process(target=do_something, args=(n,x,queue)) for n,x in enumerate(data)]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

    unsorted_result = [queue.get() for _ in processes]

    result = [i[1] for i in sorted(unsorted_result)]
    print(result)


if __name__ == '__main__':
    main()
