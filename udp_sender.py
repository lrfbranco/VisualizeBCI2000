# udp_sender.py

import socket
import time
import random

def main():
    host = 'localhost'
    port = 9999
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    # 1) Send channel properties once:
    ch_names = ['ch1','ch2','ch3']
    props = 'PROPS:' + ','.join(ch_names)
    sock.sendto(props.encode('utf-8'), (host, port))

    # 2) Loop sending random data every 100 ms
    while True:
        values = [random.random() for _ in ch_names]
        msg = 'DATA:' + ','.join(f"{v:.4f}" for v in values)
        sock.sendto(msg.encode('utf-8'), (host, port))
        print('Data sent ', values[1])
        time.sleep(0.1)

if __name__ == '__main__':
    main()
