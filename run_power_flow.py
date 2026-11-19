"""Run the included three-bus AC load-flow case and export results."""
import csv, json
from pathlib import Path
from power_flow import solve_case

ROOT=Path(__file__).resolve().parent

def main():
    case=json.loads((ROOT/"data"/"case_3bus.json").read_text(encoding="utf-8"))
    result,base=solve_case(case)
    if not result.converged:
        raise SystemExit(f"Did not converge: mismatch={result.mismatch_norm:.3e} pu")
    out=ROOT/"results"; out.mkdir(exist_ok=True)
    with (out/"bus_results.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f); w.writerow(["bus_id","voltage_pu","angle_deg","p_injection_mw","q_injection_mvar"])
        for i,bid in enumerate(result.bus_ids):
            w.writerow([bid,f"{result.voltage_magnitudes[i]:.8f}",
                        f"{result.voltage_angles_deg[i]:.8f}",
                        f"{result.p_injections[i]*base:.8f}",f"{result.q_injections[i]*base:.8f}"])
    with (out/"branch_flows.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(result.branch_flows[0])); w.writeheader(); w.writerows(result.branch_flows)
    summary={"case_name":case["case_name"],"base_mva":base,"converged":result.converged,
             "iterations":result.iterations,"maximum_mismatch_pu":result.mismatch_norm,
             "total_active_loss_mw":result.total_loss_pu.real*base,
             "total_reactive_loss_mvar":result.total_loss_pu.imag*base}
    (out/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(f"Case: {case['case_name']}")
    print(f"Converged: {result.converged} in {result.iterations} iterations")
    print(f"Maximum mismatch: {result.mismatch_norm:.3e} pu")
    print(f"Total active losses: {summary['total_active_loss_mw']:.6f} MW")
    print(f"Total reactive losses: {summary['total_reactive_loss_mvar']:.6f} MVAr")
    print("Bus voltage results:")
    for i,bid in enumerate(result.bus_ids):
        print(f"Bus {bid}: {result.voltage_magnitudes[i]:.6f} pu, {result.voltage_angles_deg[i]:.4f} deg")

if __name__=="__main__": main()
