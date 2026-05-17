from alis_mock_generator import ALISMockConfig, generate_alis_mock_cube
from alis_mock_plot import plot_cube_single_band

cfg = ALISMockConfig(nx=1024, ny=1024)
cube, wl, ice_map, bd_map, meta, gd = generate_alis_mock_cube(config=cfg)

plot_cube_single_band(cube, target_wavelength_nm=1500)

print(cube.shape)