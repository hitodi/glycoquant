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


def test_config_chemistry_propagates_to_generation():
    """라벨을 바꾸면 생성된 후보의 '중성질량'이 실제로 달라져야 함(일반화 핵심)."""
    from glyco import compositions
    proa = cfgmod.load()
    cands_proa = compositions.from_config(proa)
    # 2AB 설정 로드
    import glob
    ab_path = [p for p in CONFIGS if "2ab" in os.path.basename(p)]
    if not ab_path:
        return
    ab = cfgmod.load(ab_path[0])
    cands_ab = compositions.from_config(ab)
    # 동일 조성의 중성질량을 비교 — 라벨이 다르므로 달라야 한다
    key = lambda c: (c["HexNAc"], c["Hex"], c["dHex"], c["Neu5Ac"], c["Neu5Gc"])
    mp = {key(x["composition"]): x["neutral"] for x in cands_proa}
    ma = {key(x["composition"]): x["neutral"] for x in cands_ab}
    common = set(mp) & set(ma)
    assert common, "공통 조성이 없음"
    # ProA(C13H21N3O) vs 2AB(C7H8N2O) 라벨질량 차이가 반영돼야 함
    diffs = [abs(mp[k] - ma[k]) for k in common]
    assert min(diffs) > 1.0, "라벨을 바꿔도 중성질량이 동일 — 설정이 질량엔진에 연결 안 됨"
    # 또한 2AB 후보 질량은 Chemistry(2AB) 와 일치해야 함
    chem_ab = Chemistry(ab)
    x = cands_ab[0]
    assert abs(x["neutral"] - chem_ab.neutral_mass(x["composition"])) < 1e-6


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
