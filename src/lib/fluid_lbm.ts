export interface FluidLBM {
    nx: number;
    ny: number;
    N: number;
    tau: number; // relaxation time
    omega: number; // collision frequency

    f: Float32Array;
    f_new: Float32Array;
    rho: Float32Array;
    u: Float32Array;
    v: Float32Array;
    solid_mask: Float32Array; // 0 for solid, 1 for fluid
}

// D2Q9 Lattice
const cx = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const cy = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];
const opposite = [0, 3, 4, 1, 2, 7, 8, 5, 6];
const flipped_y = [0, 1, 4, 3, 2, 8, 7, 6, 5];

export function initLBM(nx: number, ny: number, tau: number): FluidLBM {
    const N = nx * ny;
    const lbm: FluidLBM = {
        nx, ny, N, tau,
        omega: 1.0 / tau,
        f: new Float32Array(N * 9),
        f_new: new Float32Array(N * 9),
        rho: new Float32Array(N).fill(1.0),
        u: new Float32Array(N),
        v: new Float32Array(N),
        solid_mask: new Float32Array(N).fill(1)
    };

    // Initialize distributions to equilibrium
    for (let i = 0; i < N; i++) {
        for (let k = 0; k < 9; k++) {
            lbm.f[i * 9 + k] = w[k];
            lbm.f_new[i * 9 + k] = w[k];
        }
    }

    return lbm;
}

export function idx(lbm: FluidLBM, i: number, j: number): number {
    return i * lbm.ny + j;
}

export function step(lbm: FluidLBM, inletVelocity: number, bc: "free-slip" | "no-slip") {
    const { nx, ny, f, f_new, rho, u, v, solid_mask, omega } = lbm;

    // Collision & Macroscopic Update
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            const n = idx(lbm, i, j);

            if (solid_mask[n] === 0) continue;

            let local_rho = 0;
            let local_u = 0;
            let local_v = 0;

            for (let k = 0; k < 9; k++) {
                const f_val = f[n * 9 + k];
                local_rho += f_val;
                local_u += f_val * cx[k];
                local_v += f_val * cy[k];
            }

            // Simple inlet boundary (Left wall)
            if (i === 1 && solid_mask[n] !== 0) {
                local_u = inletVelocity;
                local_v = 0;
                local_rho = 1.0;
            } else {
                local_u /= local_rho;
                local_v /= local_rho;
            }

            rho[n] = local_rho;
            u[n] = local_u;
            v[n] = local_v;

            const u_sq = local_u**2 + local_v**2;

            // BGK Collision
            for (let k = 0; k < 9; k++) {
                const cu = cx[k] * local_u + cy[k] * local_v;
                const feq = w[k] * local_rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq);
                
                f[n * 9 + k] += omega * (feq - f[n * 9 + k]);

                // Force equilibrium at inlet
                if (i === 1 && solid_mask[n] !== 0) {
                    f[n * 9 + k] = feq;
                }
            }
        }
    }

    // Streaming
    for (let j = 0; j < ny; j++) {
        for (let i = 0; i < nx; i++) {
            const n = idx(lbm, i, j);

            if (solid_mask[n] === 0) continue;

            for (let k = 0; k < 9; k++) {
                let ni = i + cx[k];
                let nj = j + cy[k];

                // Zero-gradient (free-slip) neumann BC
                if (ni >= nx - 1) ni = nx - 2; // right wall
                
                if (bc === "free-slip") {
                    if (nj < 0) { // bottom wall
                        nj = 0;
                        f_new[idx(lbm, ni, nj) * 9 + flipped_y[k]] = f[n * 9 + k];
                    } else if (nj >= ny) { // top wall
                        nj = ny - 1;
                        f_new[idx(lbm, ni, nj) * 9 + flipped_y[k]] = f[n * 9 + k];
                    }
                } else if (bc === "no-slip") {
                    if (nj < 0 || nj >= ny) { // top/bottom walls
                        f_new[n * 9 + opposite[k]] = f[n * 9 + k];
                        continue;
                    }
                }

                // Bounce-back on internal solids
                if (solid_mask[idx(lbm, ni, nj)] === 0) {
                    f_new[n * 9 + opposite[k]] = f[n * 9 + k];
                } else {
                    f_new[idx(lbm, ni, nj) * 9 + k] = f[n * 9 + k];
                }
            }
        }
    }

    // Swap buffers
    lbm.f = f_new;
    lbm.f_new = f;
}