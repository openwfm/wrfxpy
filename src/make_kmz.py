from forecast import make_kmz
import logging
import sys

if __name__ == '__main__':

    # configure the basic logger
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    args = ' '.join(sys.argv[1:])
    make_kmz(args)

