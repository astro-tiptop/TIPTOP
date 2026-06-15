#!/usr/bin/env python3
import os
import time
import argparse
import statistics
import json
import platform
from datetime import datetime
from pathlib import Path
import numpy as np

# TIPTOP Explicit Imports
from tiptop.baseSimulation import baseSimulation
from tiptop.tiptopUtils import cpuArray
from tiptop._version import __version__

try:
    import cupy as cp
except ImportError:
    cp = None

def get_system_info():
    """Extracts hardware and OS information to contextualize the benchmark."""
    info = {
        "hostname": platform.node(),
        "os": f"{platform.system()} {platform.release()}",
        "cpu_cores": os.cpu_count(),
        "gpu_enabled": False,
        "gpu_name": "None"
    }
    
    # Try to extract GPU hardware info if CuPy is available and active
    if cp is not None:
        try:
            device_id = cp.cuda.Device().id
            props = cp.cuda.runtime.getDeviceProperties(device_id)
            info["gpu_enabled"] = True
            info["gpu_name"] = props['name'].decode('utf-8')
        except Exception:
            pass
            
    return info

def _synchronize_gpu():
    """Ensures accurate timing by waiting for the GPU to finish its pending tasks."""
    if cp is not None:
        try:
            cp.cuda.runtime.deviceSynchronize()
        except Exception:
            pass

def _timed_call(func):
    """Executes a function and safely times its execution."""
    _synchronize_gpu()
    t0 = time.perf_counter()
    func()
    _synchronize_gpu()
    return time.perf_counter() - t0

def run_simulation_case(config_dir: Path, case_name: str, output_dir: Path):
    """Initializes and runs a single TIPTOP simulation."""
    sim = baseSimulation(
        path=str(config_dir),
        parametersFile=case_name,
        outputDir=str(output_dir),
        outputFile=f"bench_{case_name}_out",
        doConvolve=True,
        doPlot=False,
        addSrAndFwhm=True,
        verbose=False
    )
    # We pass astIndex=None to run the full field evaluation
    sim.doOverallSimulation(astIndex=None)
    sim.computeMetrics()
    return sim

def main():
    parser = argparse.ArgumentParser(description="TIPTOP Performance and Regression Benchmark")
    parser.add_argument("--config-dir", type=str, default="tiptop/perfTest", help="Directory containing .ini files")
    parser.add_argument("--output-dir", type=str, default="/tmp", help="Directory for temporary outputs")
    parser.add_argument("--repeats", type=int, default=3, help="Number of timed runs per INI")
    parser.add_argument("--warmups", type=int, default=1, help="Number of warm-up runs (not timed)")
    parser.add_argument("--save-json", type=str, default="benchmark_history.json", help="Path to append JSON results")
    args = parser.parse_args()

    config_dir = Path(args.config_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # System Info Extraction
    sys_info = get_system_info()

    # Deliberately skipped files
    skip_list = ["MAVIStest.ini"]

    if not config_dir.exists():
        print(f"ERROR: The directory {config_dir} does not exist.")
        return

    # Find all .ini files in the directory
    ini_files = sorted([f for f in config_dir.glob("*.ini") if f.name not in skip_list])
    
    if not ini_files:
        print(f"No .ini files found in {config_dir}")
        return

    # Print Header with Hardware Info
    print("="*60)
    print(" TIPTOP BENCHMARK RUNNER")
    print("="*60)
    print(f" Hostname   : {sys_info['hostname']}")
    print(f" OS         : {sys_info['os']} ({sys_info['cpu_cores']} cores)")
    if sys_info['gpu_enabled']:
        print(f" GPU        : \033[92mENABLED\033[0m - {sys_info['gpu_name']}")
    else:
        print(f" GPU        : \033[91mDISABLED\033[0m (Running on CPU)")
    print("-" * 60)
    print(f" Configurations : {len(ini_files)} files in {config_dir.name}/")
    print(f" Execution      : {args.warmups} warmups, {args.repeats} repeats")
    print("="*60 + "\n")

    results = []

    for ini_path in ini_files:
        case_name = ini_path.stem
        print(f" Benchmarking: {case_name:<20} ... ", end="", flush=True)

        try:
            # WARMUP PHASE
            for _ in range(args.warmups):
                run_simulation_case(config_dir, case_name, output_dir)

            # TIMED RUNS PHASE
            samples = []
            last_sim = None
            for _ in range(args.repeats):
                t0 = time.perf_counter()
                _synchronize_gpu()
                last_sim = run_simulation_case(config_dir, case_name, output_dir)
                _synchronize_gpu()
                dt = time.perf_counter() - t0
                samples.append(dt)

            avg_time = statistics.mean(samples)
            std_time = statistics.stdev(samples) if len(samples) > 1 else 0.0

            # Safe extraction of physical metrics (taking the first science target)
            sr_val = float(np.asarray(cpuArray(last_sim.sr)).ravel()[0]) if last_sim.sr else 0.0
            fwhm_val = float(np.asarray(cpuArray(last_sim.fwhm)).ravel()[0]) if last_sim.fwhm else 0.0

            results.append({
                "case": case_name,
                "status": "OK",
                "avg_time_s": round(avg_time, 4),
                "std_time_s": round(std_time, 4),
                "sr": round(sr_val, 5),
                "fwhm_mas": round(fwhm_val, 3)
            })
            print(f"OK ({avg_time:.2f}s)")

        except Exception as e:
            results.append({
                "case": case_name,
                "status": "FAILED",
                "avg_time_s": 0.0,
                "std_time_s": 0.0,
                "sr": 0.0,
                "fwhm_mas": 0.0,
                "error": str(e).split('\n')[0] # Capture only the first line of the exception
            })
            print("FAILED")

    # ==========================================
    # FINAL ASCII TABLE OUTPUT
    # ==========================================
    print("\n\n" + "="*90)
    print(f"{'INSTRUMENT (.ini)':<20} | {'STATUS':<8} | {'AVG TIME':<15} | {'STREHL RATIO':<12} | {'FWHM (mas)':<10}")
    print("-" * 90)

    for r in results:
        if r["status"] == "OK":
            time_str = f"{r['avg_time_s']:.2f}s ± {r['std_time_s']:.2f}s"
            sr_str = f"{r['sr']:.4f}"
            fwhm_str = f"{r['fwhm_mas']:.2f}"
            print(f"{r['case']:<20} | \033[92m{r['status']:<8}\033[0m | {time_str:<15} | {sr_str:<12} | {fwhm_str:<10}")
        else:
            print(f"{r['case']:<20} | \033[91m{r['status']:<8}\033[0m | {r['error'][:40]:<45}")
            
    print("=" * 90 + "\n")

    # ==========================================
    # JSON EXPORT
    # ==========================================
    if args.save_json:
        export_data = {
            "timestamp": datetime.now().isoformat(),
            "tiptop_version": __version__,
            "system_info": sys_info,
            "benchmark_config": {
                "warmups": args.warmups,
                "repeats": args.repeats
            },
            "results": results
        }
        
        # Append to existing file if it exists, otherwise create a new list
        history = []
        if os.path.exists(args.save_json):
            try:
                with open(args.save_json, "r") as f:
                    history = json.load(f)
            except json.JSONDecodeError:
                pass # File is corrupted or empty, start fresh
                
        history.append(export_data)
        
        with open(args.save_json, "w") as f:
            json.dump(history, f, indent=4)
        print(f"Benchmark results successfully appended to {args.save_json}")

if __name__ == "__main__":
    main()
