"""Newton-Raphson AC power-flow solver for slack, PV and PQ buses."""
from dataclasses import dataclass
from typing import Any
import numpy as np

@dataclass(frozen=True)
class Bus:
    bus_id: int
    bus_type: str
    p_spec: float
    q_spec: float
    vm_spec: float
    va_spec_deg: float = 0.0

@dataclass(frozen=True)
class Branch:
    from_bus: int
    to_bus: int
    r_pu: float
    x_pu: float
    b_pu: float = 0.0

@dataclass
class PowerFlowResult:
    converged: bool
    iterations: int
    mismatch_norm: float
    bus_ids: list
    voltage_magnitudes: np.ndarray
    voltage_angles_deg: np.ndarray
    p_injections: np.ndarray
    q_injections: np.ndarray
    branch_flows: list
    total_loss_pu: complex
    y_bus: np.ndarray

class PowerFlowError(ValueError):
    pass

def parse_case(case: dict[str, Any]):
    try:
        base = float(case["base_mva"])
        buses = [Bus(int(x["bus_id"]), str(x["type"]).upper(),
                     float(x.get("p_pu", 0)), float(x.get("q_pu", 0)),
                     float(x.get("vm_pu", 1)), float(x.get("va_deg", 0)))
                 for x in case["buses"]]
        branches = [Branch(int(x["from_bus"]), int(x["to_bus"]),
                           float(x["r_pu"]), float(x["x_pu"]),
                           float(x.get("b_pu", 0))) for x in case["branches"]]
    except (KeyError, TypeError, ValueError) as exc:
        raise PowerFlowError(f"Invalid case data: {exc}") from exc
    if not np.isfinite(base) or base <= 0 or not buses:
        raise PowerFlowError("A positive base_mva and at least one bus are required.")
    ids = [x.bus_id for x in buses]
    if len(ids) != len(set(ids)):
        raise PowerFlowError("Bus IDs must be unique.")
    if sum(x.bus_type == "SLACK" for x in buses) != 1:
        raise PowerFlowError("Exactly one slack bus is required.")
    if any(x.bus_type not in {"SLACK", "PV", "PQ"} for x in buses):
        raise PowerFlowError("Bus types must be SLACK, PV or PQ.")
    if any(not np.all(np.isfinite([x.p_spec, x.q_spec, x.vm_spec, x.va_spec_deg]))
           or x.vm_spec <= 0 for x in buses):
        raise PowerFlowError("Bus values must be finite and voltage setpoints positive.")
    known = set(ids)
    for x in branches:
        if x.from_bus not in known or x.to_bus not in known or x.from_bus == x.to_bus:
            raise PowerFlowError("Invalid branch bus reference.")
        if not np.all(np.isfinite([x.r_pu, x.x_pu, x.b_pu])) or abs(complex(x.r_pu,x.x_pu)) == 0:
            raise PowerFlowError("Invalid branch impedance or charging.")
    return buses, branches, base

def build_ybus(buses, branches):
    idx = {x.bus_id:i for i,x in enumerate(buses)}
    y = np.zeros((len(buses),len(buses)), dtype=complex)
    for line in branches:
        i,j = idx[line.from_bus],idx[line.to_bus]
        ys = 1/complex(line.r_pu,line.x_pu)
        yh = 1j*line.b_pu/2
        y[i,i] += ys+yh; y[j,j] += ys+yh
        y[i,j] -= ys; y[j,i] -= ys
    return y

def calculate_power(voltage, ybus):
    s = voltage*np.conj(ybus@voltage)
    return s.real,s.imag

def _jacobian(vm,va,p,q,g,b,non_slack,pq):
    na,nv = len(non_slack),len(pq)
    jmat = np.zeros((na+nv,na+nv))
    for row,i in enumerate(non_slack):
        for col,k in enumerate(non_slack):
            if i == k: jmat[row,col] = -q[i]-b[i,i]*vm[i]**2
            else:
                d=va[i]-va[k]
                jmat[row,col]=vm[i]*vm[k]*(g[i,k]*np.sin(d)-b[i,k]*np.cos(d))
        for col,k in enumerate(pq):
            if i == k: jmat[row,na+col]=p[i]/vm[i]+g[i,i]*vm[i]
            else:
                d=va[i]-va[k]
                jmat[row,na+col]=vm[i]*(g[i,k]*np.cos(d)+b[i,k]*np.sin(d))
    for row,i in enumerate(pq):
        for col,k in enumerate(non_slack):
            if i == k: jmat[na+row,col]=p[i]-g[i,i]*vm[i]**2
            else:
                d=va[i]-va[k]
                jmat[na+row,col]=-vm[i]*vm[k]*(g[i,k]*np.cos(d)+b[i,k]*np.sin(d))
        for col,k in enumerate(pq):
            if i == k: jmat[na+row,na+col]=q[i]/vm[i]-b[i,i]*vm[i]
            else:
                d=va[i]-va[k]
                jmat[na+row,na+col]=vm[i]*(g[i,k]*np.sin(d)-b[i,k]*np.cos(d))
    return jmat

def solve_power_flow(buses,branches,tolerance=1e-10,max_iterations=30):
    if tolerance<=0 or max_iterations<1: raise PowerFlowError("Invalid solver settings.")
    if len({x.bus_id for x in buses}) != len(buses) or sum(x.bus_type=="SLACK" for x in buses)!=1:
        raise PowerFlowError("Require unique bus IDs and exactly one slack bus.")
    ybus=build_ybus(buses,branches); g,b=ybus.real,ybus.imag; n=len(buses)
    slack=next(i for i,x in enumerate(buses) if x.bus_type=="SLACK")
    non_slack=np.array([i for i in range(n) if i!=slack],dtype=int)
    pq=np.array([i for i,x in enumerate(buses) if x.bus_type=="PQ"],dtype=int)
    ps=np.array([x.p_spec for x in buses]); qs=np.array([x.q_spec for x in buses])
    vm=np.array([x.vm_spec for x in buses],dtype=float)
    va=np.radians(np.array([x.va_spec_deg for x in buses],dtype=float))
    norm=float("inf"); converged=False; iterations=0
    for it in range(max_iterations+1):
        v=vm*np.exp(1j*va); pc,qc=calculate_power(v,ybus)
        mismatch=np.r_[ps[non_slack]-pc[non_slack],qs[pq]-qc[pq]]
        norm=float(np.max(np.abs(mismatch))) if mismatch.size else 0.0
        if norm<tolerance: converged=True; iterations=it; break
        if it==max_iterations: iterations=it; break
        try: dx=np.linalg.solve(_jacobian(vm,va,pc,qc,g,b,non_slack,pq),mismatch)
        except np.linalg.LinAlgError as exc: raise PowerFlowError("Singular Jacobian.") from exc
        va[non_slack]+=dx[:len(non_slack)]
        if len(pq):
            vm[pq]+=dx[len(non_slack):]
            if np.any(vm[pq]<=0): raise PowerFlowError("Non-positive voltage update.")
    v=vm*np.exp(1j*va); pc,qc=calculate_power(v,ybus)
    idx={x.bus_id:i for i,x in enumerate(buses)}
    flows=[]; loss_total=0j
    for line in branches:
        i,k=idx[line.from_bus],idx[line.to_bus]
        ys=1/complex(line.r_pu,line.x_pu); yh=1j*line.b_pu/2
        sf=v[i]*np.conj((v[i]-v[k])*ys+v[i]*yh)
        st=v[k]*np.conj((v[k]-v[i])*ys+v[k]*yh)
        loss=sf+st; loss_total+=loss
        flows.append({"from_bus":line.from_bus,"to_bus":line.to_bus,
                      "p_from_pu":float(sf.real),"q_from_pu":float(sf.imag),
                      "p_to_pu":float(st.real),"q_to_pu":float(st.imag),
                      "p_loss_pu":float(loss.real),"q_loss_pu":float(loss.imag)})
    return PowerFlowResult(converged,iterations,norm,[x.bus_id for x in buses],
                           vm.copy(),np.degrees(va),pc,qc,flows,loss_total,ybus)

def solve_case(case,tolerance=1e-10,max_iterations=30):
    buses,branches,base=parse_case(case)
    return solve_power_flow(buses,branches,tolerance,max_iterations),base
