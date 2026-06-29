"""설정 로딩 + 설정별 Chemistry 가 진단이온을 재현하는지."""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from glyco import config as cfgmod
from glyco.chem import Chemistry

CONFIGS = glob.glob(os.path.join(os.path.dirname(cfgmod.CONFIG_DIR), "configs", "*.yaml"))


def test_all_configs_load():
    assert CONFIGS, "configs/*.yaml 가 없습니다"
    for path in CONFIGS:
        cfg = cfgmod.load(path)
        chem = Chemistry(cfg)
        # 모든 진단이온이 양수 m/z 로 계산되는지
        for name, mz in chem.diagnostic_table():
            assert mz > 0, f"{path}: {name} m/z={mz}"


def test_default_proa_diagnostics():
    cfg = cfgmod.load()           # 기본 ProA
    chem = Chemistry(cfg)
    table = dict(chem.diagnostic_table())
    assert abs(table["HexNAc"] - 204.0867) < 0.01
    assert abs(table["ProA-HexNAc"] - 441.2708) < 0.01


if __name__ == "__main__":
    test_all_configs_load()
    test_default_proa_diagnostics()
    print("설정 테스트 통과 ✓")
