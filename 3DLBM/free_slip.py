import numpy as np
import os

# D3Q19 Lattice Constants
cx = np.array([0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0])
cy = np.array([0, 0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1])
cz = np.array([0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1])

w = np.array([1/3] + [1/18]*6 + [1/36]*12)

# Reversal Arrays
opposite = np.array([0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17])
# Specular reflection for Free-Slip (Inverts Y-velocity only)
flipped_y = np.array([0, 1, 2, 4, 3, 5, 6, 9, 10, 7, 8, 11, 12, 13, 14, 18, 17, 16, 15])

class FluidLBM:
    def __init__(self, nx, ny, nz, tau):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.tau = tau
        self.omega = 1.0 / tau
        
        self.f = np.zeros((19, nx, ny, nz))
        self.f_new = np.zeros((19, nx, ny, nz))
        self.rho = np.ones((nx, ny, nz))
        self.u = np.zeros((nx, ny, nz))
        self.v = np.zeros((nx, ny, nz))
        self.w_vel = np.zeros((nx, ny, nz))
        self.solid_mask = np.ones((nx, ny, nz), dtype=int) 

def initLBM(nx, ny, nz, tau, u_inlet):
    lbm = FluidLBM(nx, ny, nz, tau)
    # Initialize entire domain with freestream velocity
    u_sq = u_inlet**2
    for k in range(19):
        cu = cx[k] * u_inlet
        f_eq = w[k] * 1.0 * (1.0 + 3.0*cu + 4.5*cu**2 - 1.5*u_sq)
        lbm.f[k, :, :, :] = f_eq
    return lbm

def step(lbm, inletVelocity):
    # Collision & Macroscopic Update
    lbm.rho = np.sum(lbm.f, axis=0)
    lbm.u = np.sum(lbm.f * cx[:, None, None, None], axis=0) / lbm.rho
    lbm.v = np.sum(lbm.f * cy[:, None, None, None], axis=0) / lbm.rho
    lbm.w_vel = np.sum(lbm.f * cz[:, None, None, None], axis=0) / lbm.rho

    lbm.u[lbm.solid_mask == 0] = 0
    lbm.v[lbm.solid_mask == 0] = 0
    lbm.w_vel[lbm.solid_mask == 0] = 0

    u_sq = lbm.u**2 + lbm.v**2 + lbm.w_vel**2
    f_eq = np.zeros_like(lbm.f)
    
    # BGK Collision
    for k in range(19):
        cu = cx[k]*lbm.u + cy[k]*lbm.v + cz[k]*lbm.w_vel
        f_eq[k] = w[k] * lbm.rho * (1.0 + 3.0*cu + 4.5*cu**2 - 1.5*u_sq)
        
    lbm.f_new = lbm.f + lbm.omega * (f_eq - lbm.f)

    # Force equilibrium at inlet
    u_inlet_sq = inletVelocity**2
    for k in range(19):
        cu_inlet = cx[k]*inletVelocity
        f_inlet_eq = w[k] * 1.0 * (1.0 + 3.0*cu_inlet + 4.5*cu_inlet**2 - 1.5*u_inlet_sq)
        lbm.f_new[k, 0, :, :] = f_inlet_eq

    # Bounce-back on internal solids
    for k in range(19):
        lbm.f_new[k, lbm.solid_mask == 0] = lbm.f[opposite[k], lbm.solid_mask == 0]

    # Streaming
    for k in range(19):
        lbm.f[k] = np.roll(lbm.f_new[k], shift=(cx[k], cy[k], cz[k]), axis=(0, 1, 2))

    # Domain Boundary Conditions 
    # Y-axis: Free-Slip (Specular Reflection using flipped_y)
    # Z-axis: Periodic (np.roll inherently handles this)
    for k in range(19):
        if cy[k] > 0: lbm.f[k, :, 0, :] = lbm.f_new[flipped_y[k], :, 0, :]
        if cy[k] < 0: lbm.f[k, :, -1, :] = lbm.f_new[flipped_y[k], :, -1, :]

    # Zero-gradient convective outlet
    lbm.f[:, -1, :, :] = lbm.f[:, -2, :, :]

def exportToVTK(domain, step_num, out_dir="vtk_output"):
    os.makedirs(out_dir, exist_ok=True)
    filename = os.path.join(out_dir, f"flow_{step_num:05d}.vtk")
    
    with open(filename, 'w') as f_out:
        f_out.write("# vtk DataFile Version 3.0\n")
        f_out.write("LBM D3Q19 Velocity Field\nASCII\nDATASET STRUCTURED_POINTS\n")
        f_out.write(f"DIMENSIONS {domain.nx} {domain.ny} {domain.nz}\n")
        f_out.write("ORIGIN 0 0 0\nSPACING 1 1 1\n")
        f_out.write(f"POINT_DATA {domain.nx * domain.ny * domain.nz}\n")
        f_out.write("VECTORS velocity float\n")
        
        vel = np.column_stack((
            domain.u.flatten('F'), 
            domain.v.flatten('F'), 
            domain.w_vel.flatten('F')
        ))
        np.savetxt(f_out, vel, fmt='%.5f')

if __name__ == "__main__":
    # Parameters scaled for Re = 400 stability
    nx, ny, nz = 300, 100, 30
    lbm_tau = 0.5225 
    lbm_inlet_velocity = 0.1
    niters = 25000
    save_freq = 100
    
    obstacle_size = 15
    cx_pos, cy_pos = 60, 50.1 
    
    print(f"Initializing {nx}x{ny}x{nz} domain...")
    lbm_domain = initLBM(nx, ny, nz, lbm_tau, lbm_inlet_velocity)
    
    x, y, z = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing='ij')
    cylinder_mask = (x - cx_pos)**2 + (y - cy_pos)**2 <= obstacle_size**2
    lbm_domain.solid_mask[cylinder_mask] = 0
    
    print("Starting simulation... (Press Ctrl+C to abort)")
    try:
        for current_frame in range(niters):
            step(lbm_domain, lbm_inlet_velocity)
            
            if current_frame % save_freq == 0:
                max_speed = np.max(np.sqrt(lbm_domain.u**2 + lbm_domain.v**2 + lbm_domain.w_vel**2))
                print(f"Frame {current_frame:05d}/{niters} | Max Velocity: {max_speed:.4f}")
                exportToVTK(lbm_domain, current_frame)
                
    except KeyboardInterrupt:
        print("\nSimulation interrupted by user. Saving final frame...")
        exportToVTK(lbm_domain, current_frame)
        
    print("Done. You can now load the vtk_output folder in ParaView.")
