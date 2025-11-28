import os
import subprocess
import time
import shutil

def benchmark():
    print("Compiling...")
    try:
        subprocess.run(["ninja", "-C", "build"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print("Compilation failed!")
        return

    # Ensure no stale restart files pollute timing
    for root_file in ["xtbrestart"]:
        if os.path.exists(root_file):
            os.remove(root_file)

    xtb_bin = os.path.abspath("build/xtb")
    input_file = os.path.abspath("9NY3.pdb")
    work_dir = "bench_9ny3_tmp"

    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)

    os.chdir(work_dir)
    shutil.copy(input_file, "9NY3.pdb")

    times = []
    print("Running benchmark on 9NY3.pdb...")
    for i in range(3):
        # Clean up potential restart files
        for f in ["xtbrestart", "xtbtopo.mol", "charges", "wbo", "xtbopt.log", "xtbopt.coord", "xtbla.restart"]:
            if os.path.exists(f):
                os.remove(f)

        start = time.time()
        try:
            # Added --uhf 1 for RCSB 9NY3.pdb compatibility
            subprocess.run([xtb_bin, "9NY3.pdb", "--gfn", "2", "--uhf", "1"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Run {i+1} failed!")
            break
        end = time.time()
        duration = end - start
        print(f"Run {i+1}: {duration:.4f} s")
        times.append(duration)

    if times:
        avg = sum(times) / len(times)
        print(f"Average time: {avg:.4f} s")

    os.chdir("..")
    shutil.rmtree(work_dir)

if __name__ == "__main__":
    benchmark()
