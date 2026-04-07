// Macroscopic eulerian FDM, Chorin's Projection Method, MAC staggered grid
// incompressible, inviscid, constant density fluid, neumann BC

// todo: refactor monolithic code into separate files

export interface Fluid {
    nx: number;
    ny: number;
    h: number;
    N: number;

    u: Float32Array;
    v: Float32Array;
    pressure: Float32Array;
    solid_mask: Float32Array;
    fluid_density: Float32Array;
}

export function initFluid(nx: number, ny: number, h: number): Fluid {
    const N = (nx + 2) * (ny + 2);
    const fluid: Fluid = {
        // density: density, // assuming constant density throughout
        nx: nx + 2, // adding buffer cells
        ny: ny + 2, // for boundary handling
        h: h, // cell size
        N: N, // total cells
        u: new Float32Array(N), // horizontal velocity
        v: new Float32Array(N), // vertical velocity
        pressure: new Float32Array(N),
        solid_mask: new Float32Array(N).fill(1), // 0 for solid, 1 for fluid
        fluid_density: new Float32Array(N).fill(1) // clear fluid
    };
    return fluid;
}

export function idx(fluid: Fluid, i: number, j: number): number {
    return i * fluid.ny + j;
}

// Euler ODE approximation (y_1 = y_0 + f(y_0) * (t_1 - t_0))
// export function integrateGravity(fluid: Fluid, g: number, dt: number) {
//     // todo: add other external forces

//     for (let i = 1; i < fluid.nx - 1; i++) {
//         for (let j = 1; j < fluid.ny - 1; j++) { // leave out buffer cells
//             if (fluid.solid_mask[idx(fluid, i, j)] !== 0 && fluid.solid_mask[idx(fluid, i, j-1)]) { // skip if cell b'ary above and below is solid
//                 fluid.v[idx(fluid, i, j)] += g * dt; // apply gravity to vertical velocity
//             }
//         }
//     }
// }

// Viscous diffusion (u += \nu dt \del^2u)
export function diffuse(fluid: Fluid, nu: number, dt: number, iters: number) {
    if (nu === 0) return;   // skip if inviscid

    const coeff = nu * dt / (fluid.h * fluid.h);

    const u0 = new Float32Array(fluid.u);
    const v0 = new Float32Array(fluid.v);

    for (let iter = 0; iter < iters; iter++) {
        for (let i = 1; i < fluid.nx - 1; i++) {
            for (let j = 1; j < fluid.ny - 1; j++) { // leave out buffer cells
                
                // Diffuse horizontal velocity (u)
                if (fluid.solid_mask[idx(fluid, i, j)] !== 0 && fluid.solid_mask[idx(fluid, i - 1, j)] !== 0) {
                    fluid.u[idx(fluid, i, j)] = (u0[idx(fluid, i, j)] + coeff * (
                        fluid.u[idx(fluid, i - 1, j)] + fluid.u[idx(fluid, i + 1, j)] +
                        fluid.u[idx(fluid, i, j - 1)] + fluid.u[idx(fluid, i, j + 1)]
                    )) / (1 + 4 * coeff);
                }
                
                // Diffuse vertical velocity (v)
                if (fluid.solid_mask[idx(fluid, i, j)] !== 0 && fluid.solid_mask[idx(fluid, i, j - 1)] !== 0) {
                    fluid.v[idx(fluid, i, j)] = (v0[idx(fluid, i, j)] + coeff * (
                        fluid.v[idx(fluid, i - 1, j)] + fluid.v[idx(fluid, i + 1, j)] +
                        fluid.v[idx(fluid, i, j - 1)] + fluid.v[idx(fluid, i, j + 1)]
                    )) / (1 + 4 * coeff);
                }
            }
        }
    }

}

// Incompressibility divergence-free condition (gauss-siedel SOR)
export function solvePressure(fluid: Fluid, iters: number, dt: number, over_relaxation: number) {
    for (let iter = 0; iter < iters; iter++) {
        for (let i = 1; i < fluid.nx - 1; i++) {
            for (let j = 1; j < fluid.ny - 1; j++) { // leave out buffer cells
                if (fluid.solid_mask[idx(fluid, i, j)] !== 0) { // only solve for fluid cells
                    
                    const div = (fluid.u[idx(fluid, i+1, j)] - fluid.u[idx(fluid, i, j)]) + // raw difference instead of proper divergenece
                                (fluid.v[idx(fluid, i, j+1)] - fluid.v[idx(fluid, i, j)]);  // which would require division by h

                    const sx_left   = fluid.solid_mask[idx(fluid, i-1, j)]; // check if
                    const sx_right  = fluid.solid_mask[idx(fluid, i+1, j)]; // neighbouring
                    const sy_down   = fluid.solid_mask[idx(fluid, i, j-1)]; // cells are
                    const sy_up     = fluid.solid_mask[idx(fluid, i, j+1)]; // solid or fluid
                    
                    const s_coeff = sx_left + sx_right + sy_down + sy_up; // fluid neighbours (out of 4)
                    if (s_coeff === 0) continue; // skip if all neighbors are solid

                    let p = (-div / s_coeff) * over_relaxation; // average pressure contrib from neighbour cells + SOR
                    // fluid.pressure[idx(fluid, i, j)] += p * fluid.density * fluid.h / dt; // convert to pressure units & accumulate over iters
                    fluid.pressure[idx(fluid, i, j)] += p * fluid.h / dt;

                    if (sx_left !== 0) fluid.u[idx(fluid, i, j)] -= p;     // velocity
                    if (sx_right !== 0) fluid.u[idx(fluid, i+1, j)] += p;  // correction
                    if (sy_down !== 0) fluid.v[idx(fluid, i, j)] -= p;     // from
                    if (sy_up !== 0) fluid.v[idx(fluid, i, j+1)] += p;     // pressure gradient
                }
            }
        }
    }
}

// free-slip (neumann, (u_i+1 - u_i)/h = 0) b'ary conditions (zero shear stress at b'ary)
export function applyBoundaryConditions(fluid: Fluid, bc: "free-slip" | "no-slip") {
    
    if (bc === "free-slip") { // neumann
        for (let i = 0; i < fluid.nx; i++) {
            fluid.u[idx(fluid, i, 0)] = fluid.u[idx(fluid, i, 1)]; // bottom b'ary
            fluid.u[idx(fluid, i, fluid.ny - 1)] = fluid.u[idx(fluid, i, fluid.ny - 2)]; // top b'ary
        }
    } else if (bc === "no-slip") { // dirichlet
        for (let i = 0; i < fluid.nx; i++) {
            fluid.u[idx(fluid, i, 0)] = 0; // bottom b'ary
            fluid.u[idx(fluid, i, fluid.ny - 1)] = 0; // top b'ary
        }
    }
    
    for (let j = 0; j < fluid.ny; j++) {
        fluid.v[idx(fluid, 0, j)] = fluid.v[idx(fluid, 1, j)]; // left b'ary
        fluid.v[idx(fluid, fluid.nx - 1, j)] = fluid.v[idx(fluid, fluid.nx - 2, j)]; // right b'ary
        fluid.u[idx(fluid, fluid.nx - 1, j)] = fluid.u[idx(fluid, fluid.nx - 2, j)];
    }
}

// parameter value at arbitrary position
export function arbitraryPosnParameter(fluid: Fluid, x: number, y: number, param: String): number {
    x = Math.max(fluid.h, Math.min(x, (fluid.nx - 2) * fluid.h)); // clamp to
    y = Math.max(fluid.h, Math.min(y, (fluid.ny - 2) * fluid.h)); // valid grid range

    let dx = 0, dy = 0, p;
    switch (param) {
        case "u":
            dy = 0.5 * fluid.h; // u is staggered in vertical face
            p = fluid.u; break;
        case "v":
            dx = 0.5 * fluid.h; // v is staggered in horizontal face
            p = fluid.v; break;
        default:
            p = fluid.fluid_density;
            dx = 0.5 * fluid.h; // center
            dy = 0.5 * fluid.h; // of cell
    }

    const x0 = Math.floor((x - dx) / fluid.h); // integer index of grid node to the left of x
    const y0 = Math.floor((y - dy) / fluid.h); // integer index of grid node below y
    
    const frac_x_from_left = (x - dx) / fluid.h - x0; // fractional distance of x from left grid node
    const frac_y_from_left = (y - dy) / fluid.h - y0; // fractional distance of y from bottom grid node
    
    const frac_x_from_right = 1 - frac_x_from_left; // fractional distance of x from right grid node
    const frac_y_from_right = 1 - frac_y_from_left; // fractional distance of y from top grid node
    
    const x1 = x0 + 1; // integer index of grid node to the right of x
    const y1 = y0 + 1; // integer index of grid node above y

    // bilinear interpolation
    return p[idx(fluid, x0, y0)] * frac_x_from_right * frac_y_from_right + // bottom left corner weighted by top right quad area
            p[idx(fluid, x1, y0)] * frac_x_from_left * frac_y_from_right +  // bottom right corner weighted by top left quad area
            p[idx(fluid, x0, y1)] * frac_x_from_right * frac_y_from_left +  // top left corner weighted by bottom right quad area
            p[idx(fluid, x1, y1)] * frac_x_from_left * frac_y_from_left;    // top right corner weighted by bottom left quad area
}

// calculate velocities at edge centers
export function averageVelocity(fluid: Fluid, i: number, j: number, direction: String): number {
    if (direction === "u") {
        return 0.25 * (fluid.u[idx(fluid, i, j)] + fluid.u[idx(fluid, i+1, j)] +     // value of u where normally
                     fluid.u[idx(fluid, i, j-1)] + fluid.u[idx(fluid, i+1, j-1)])    // v is defined (horizontal face)
    } else {
        return 0.25 * (fluid.v[idx(fluid, i-1, j+1)] + fluid.v[idx(fluid, i, j+1)] + // value of v where normally
                         fluid.v[idx(fluid, i-1, j)] + fluid.v[idx(fluid, i, j)])    // u is defined (vertical face)
    }
}

// advection (todo: review code)
export function advect(fluid: Fluid, dt: number) {
    let newU = new Float32Array(fluid.u);
    let newV = new Float32Array(fluid.v);
    let newDensity = new Float32Array(fluid.fluid_density);

    for (let i = 1; i < fluid.nx - 1; i++) {
        for (let j = 1; j < fluid.ny - 1; j++) { // leave out buffer cells
            
            if (fluid.solid_mask[idx(fluid, i, j)] !== 0 && fluid.solid_mask[idx(fluid, i-1, j)] !== 0) { // skip if cell b'ary left and right is solid
                const x = i * fluid.h;
                const y = (j + 0.5) * fluid.h; // todo: minus?
                const u_vel = fluid.u[idx(fluid, i, j)];
                const v_vel = averageVelocity(fluid, i, j, "v");
                newU[idx(fluid, i, j)] = arbitraryPosnParameter(fluid, x - u_vel * dt, y - v_vel * dt, "u"); // backtrace u array
            }
            
            if (fluid.solid_mask[idx(fluid, i, j)] !== 0 && fluid.solid_mask[idx(fluid, i, j-1)] !== 0) { // skip if cell b'ary above and below is solid
                const x = (i + 0.5) * fluid.h; // todo: minus?
                const y = j * fluid.h;
                const u_vel = averageVelocity(fluid, i, j, "u");
                const v_vel = fluid.v[idx(fluid, i, j)];
                newV[idx(fluid, i, j)] = arbitraryPosnParameter(fluid, x - u_vel * dt, y - v_vel * dt, "v"); // backtrace v array
            }

            if (fluid.solid_mask[idx(fluid, i, j)] !== 0) { // skip if cell is solid
                const u = (fluid.u[idx(fluid, i, j)] + fluid.u[idx(fluid, i+1, j)]) * 0.5; // average u velocity at cell center
                const v = (fluid.v[idx(fluid, i, j)] + fluid.v[idx(fluid, i, j+1)]) * 0.5; // average v velocity at cell center
                newDensity[idx(fluid, i, j)] = arbitraryPosnParameter(
                    fluid, (i + 0.5) * fluid.h - u * dt, (j + 0.5) * fluid.h - v * dt, "fluid_density"
                ); // backtrace from cell center
            }
        }
    }

    fluid.u = newU;
    fluid.v = newV;
    fluid.fluid_density = newDensity;
}

// simulate step
export function step(fluid: Fluid, dt: number, pressure_iters: number, over_relaxation: number, nu: number, diffuse_iters: number, bc: "free-slip" | "no-slip") {
    // integrateGravity(fluid, g, dt);
    if (nu > 0) diffuse(fluid, nu, dt, diffuse_iters);
    fluid.pressure.fill(0); // reset pressure field or contributions accumulate in solvePressure over iters
    solvePressure(fluid, pressure_iters, dt, over_relaxation);
    applyBoundaryConditions(fluid, bc); // todo: change BC
    advect(fluid, dt);
}