"""
main.py
========
全グループのモックキューブを生成してプロットを保存するエントリーポイント。
"""

from pathlib import Path

from alis_mock_generator import ALISMockConfig, DATA_DIR, generate_alis_mock_cube
from alis_mock_plot import (
    plot_calibration,
    plot_dry_spectra_all_groups,
    plot_cube_overview,
    plot_spectra_comparison,
    plot_unit_absorption_profiles,
)

if __name__ == "__main__":
    out_dir = Path(__file__).parent / "alis_mock_output_araki"
    out_dir.mkdir(exist_ok=True)

    print("Plotting dry spectra (all groups)...")
    plot_dry_spectra_all_groups(DATA_DIR, out_dir / "dry_spectra_all_groups.png")

    print("Plotting unit absorption profiles (all groups)...")
    plot_unit_absorption_profiles(DATA_DIR, out_dir / "unit_absorption_profiles.png")

    configs = [
        ALISMockConfig(mineral_type="olivine",       grain_size="coarse",
                            ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="olivine",       grain_size="fine",
                            ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="plagioclase",   grain_size="coarse",
                            ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="plagioclase",   grain_size="fine",
                            ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="clinopyroxene", grain_size="coarse",
                            ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="clinopyroxene", grain_size="fine",
                            ice_pattern="gradient", snr=100),
        ALISMockConfig(mineral_type="mixture",       grain_size="coarse",
                            ice_pattern="patches",  snr=80),
        ALISMockConfig(mineral_type="mixture",       grain_size="fine",
                            ice_pattern="patches",  snr=80),
    ]

    for cfg in configs:
        tag = f"{cfg.mineral_type}_{cfg.grain_size}"
        print(f"\n{'='*60}")
        print(f"Generating: {cfg.mineral_type} ({cfg.grain_size})")
        print(f"{'='*60}")

        cube, wl, ice_map, bd_map, meta, gd = generate_alis_mock_cube(cfg, DATA_DIR)

        print(f"  Cube shape           : {cube.shape}")
        print(f"  Wavelength range     : {wl[0]:.0f}–{wl[-1]:.0f} nm ({len(wl)} bands)")
        print(f"  Gradient (exp / Eq1) : {meta['gradient_experimental']:.2f} / "
              f"{meta['gradient_predicted_eq1']:.2f}")
        print(f"  R(1.5μm)             : {meta['reflectance_1500']:.3f}")
        print(f"  Ice range            : {ice_map.min():.2f}–{ice_map.max():.2f} wt.%")
        print(f"  BD range             : {bd_map.min():.4f}–{bd_map.max():.4f}")
        print(f"  Dry files used       : {meta['n_dry_files']} dry, {meta['n_pairs']} pairs")

        plot_calibration(ice_map, bd_map, meta, out_dir / f"{tag}_calibration.png")
        print(f"  calibration.png -> {out_dir / f'{tag}_calibration.png'}")

        plot_cube_overview(cube, ice_map, bd_map, meta, out_dir / f"{tag}_montage.png")
        print(f"  montage.png     -> {out_dir / f'{tag}_montage.png'}")

        plot_spectra_comparison(cube, ice_map, gd, meta,
                                out_dir / f"{tag}_spectra.png")
        print(f"  spectra.png     -> {out_dir / f'{tag}_spectra.png'}")

    print(f"\nAll outputs saved to: {out_dir}/")
