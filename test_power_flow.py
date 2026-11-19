import json, sys, unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from power_flow import PowerFlowError,build_ybus,calculate_power,parse_case,solve_case

class TestPowerFlow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case=json.loads((ROOT/"data"/"case_3bus.json").read_text())
        cls.buses,cls.branches,cls.base=parse_case(cls.case)
        cls.result,_=solve_case(cls.case)

    def test_ybus_symmetric(self):
        self.assertTrue(np.allclose(build_ybus(self.buses,self.branches),
                                    build_ybus(self.buses,self.branches).T,atol=1e-12))
    def test_ybus_row_sum_equals_shunts(self):
        y=build_ybus(self.buses,self.branches); expected=np.zeros(len(self.buses),complex)
        idx={x.bus_id:i for i,x in enumerate(self.buses)}
        for line in self.branches:
            expected[idx[line.from_bus]]+=1j*line.b_pu/2
            expected[idx[line.to_bus]]+=1j*line.b_pu/2
        self.assertTrue(np.allclose(y.sum(axis=1),expected,atol=1e-12))
    def test_convergence(self):
        self.assertTrue(self.result.converged)
        self.assertLess(self.result.mismatch_norm,1e-9)
    def test_slack_and_pv_voltage_setpoints(self):
        self.assertAlmostEqual(self.result.voltage_magnitudes[0],1.04,places=10)
        self.assertAlmostEqual(self.result.voltage_magnitudes[1],1.02,places=10)
        self.assertAlmostEqual(self.result.voltage_angles_deg[0],0.0,places=10)
    def test_active_power_balance(self):
        self.assertAlmostEqual(float(np.sum(self.result.p_injections)),
                               self.result.total_loss_pu.real,places=9)
    def test_branch_losses_sum_to_network_loss(self):
        self.assertAlmostEqual(sum(x["p_loss_pu"] for x in self.result.branch_flows),
                               self.result.total_loss_pu.real,places=10)
    def test_missing_slack_rejected(self):
        bad=json.loads(json.dumps(self.case)); bad["buses"][0]["type"]="PQ"
        with self.assertRaises(PowerFlowError): parse_case(bad)
    def test_power_values_finite(self):
        v=self.result.voltage_magnitudes*np.exp(1j*np.radians(self.result.voltage_angles_deg))
        p,q=calculate_power(v,self.result.y_bus)
        self.assertTrue(np.all(np.isfinite(p))); self.assertTrue(np.all(np.isfinite(q)))

if __name__=="__main__": unittest.main()
