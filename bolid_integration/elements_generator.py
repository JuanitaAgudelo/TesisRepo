import numpy as np
from tqdm import tqdm

N = int(1e8)

for start in tqdm(range(10)):
    file_elements = f"datos/random_uniform_elements/state_elements_{start}.dat"
    print('Generating orbital elements uniform distribution')
    a_uniform = np.random.uniform(0, 2, N)
    e_uniform = np.random.uniform(0, 1, N)
    i_uniform = np.random.uniform(0, 180, N)
    w_uniform = np.random.uniform(0, 360, N)
    Omega_uniform = np.random.uniform(0, 360, N)
    M_uniform = np.random.uniform(0, 360, N)

    orbital_elements = np.column_stack((a_uniform, e_uniform, i_uniform, w_uniform, Omega_uniform, M_uniform))

    np.savetxt(file_elements, orbital_elements)