from tqdm import tqdm
from optparse import OptionParser

import multiprocessing, os, signal, subprocess, time

def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)

def count_mutants(nameFile):
    output = subprocess.Popen(['wc', '-l', f'{nameFile}.seq'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = output.communicate()
    mutants = stdout.split()[0].decode('utf-8')
    return int(mutants)

def run_tango(nameFile):
    subprocess.run(f'~/bin/tango -inputfile={nameFile}.seq > log.out', shell=True)

def clean_tmp(nameFile):
    subprocess.run(f'cp *_aggregation.txt* ../../FINAL_RESULTS/{nameFile}_aggregation.txt', shell=True)
    subprocess.run('rm *.txt', shell=True)

def main(folder):
    curr_dir = os.chdir(folder)
    nameFile = folder.split('/')[-1]
    if not os.path.exists(f'../../FINAL_RESULTS/{nameFile}_aggregation.txt'):
        if count_mutants(nameFile) > 10000:
            run_tango(nameFile)
            clean_tmp(nameFile)
        else :
            print(f'{nameFile} need to be rerun on cluster')



if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", "--input_directory", dest="curr_path", default="None", help="[Required] Provide the input input_directory")
    (options, args) = parser.parse_args()
    curr_path = options.curr_path

    myfolders =  [ os.path.join(curr_path, folder) for folder in sorted(os.listdir(curr_path)) ]
    p = multiprocessing.Pool(initializer=init_worker, initargs=(tqdm.get_lock(),), processes=2)
    try:
        pbar = tqdm(myfolders, maxinterval=1.0, miniters=1, desc="Terminated Tango: ", bar_format="{desc}:{percentage:3.0f}%|{bar}|")
        for _, result in enumerate(p.imap_unordered(main, myfolders, chunksize=1)):
            pbar.update(1)  # Everytime the iteration finishes, update the global progress bar

        pbar.close()
        p.close()
        p.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt, terminating workers.")
        pbar.close()
        p.terminate()
        p.join()
        exit(1)
