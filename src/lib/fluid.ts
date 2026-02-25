// FDM + Chorin's Projection
// MAC staggered grid

// todo: refactor monolithic code into separate files

export interface Fluid {
    density: number;
    nx: number;
    ny: number;
    h: number;
    N: number;

    u: Float32Array;
    v: Float32Array;
    pressure: Float32Array;
    solid_mask: Float32Array;
    smoke_density: Float32Array;
}

export function initFluid(density: number, nx: number, ny: number, h: number): Fluid {
    const N = (nx + 2) * (ny + 2);
    const fluid: Fluid = {
        density: density, // assuming constant density throughout
        nx: nx + 2, // adding buffer cells
        ny: ny + 2, // for boundary handling
        h: h, // cell size
        N: N, // total cells
        u: new Float32Array(N), // horizontal velocity
        v: new Float32Array(N), // vertical velocity
        pressure: new Float32Array(N), // pressure field
        solid_mask: new Float32Array(N), // 0 for solid, 1 for fluid
        smoke_density: new Float32Array(N).fill(1) // clear fluid
    };
    return fluid;
}

export function idx(fluid: Fluid, i: number, j: number): number {
    return i * fluid.ny + j;
}

// Euler ODE approximation
export function integrateGravity(fluid: Fluid, g: number, dt: number) {
    for (let i = 1; i < fluid.nx - 1; i++) {
        for (let j = 1; j < fluid.ny - 1; j++) { // leave out buffer cells
            if (fluid.solid_mask[idx(fluid, i, j)] !== 0 && fluid.solid_mask[idx(fluid, i, j-1)]) { // skip if cell b'ary above and below is solid
                fluid.v[idx(fluid, i, j)] += g * dt; // apply gravity to vertical velocity
            }
        }
    }
}

// Incompressibility divergence-free condition
export function solvePressure(fluid: Fluid, iters: number, dt: number, over_relaxation: number) {
    for (let iter = 0; iter < iters; iter++) {
        for (let i = 1; i < fluid.nx - 1; i++) {
            for (let j = 1; j < fluid.ny - 1; j++) { // leave out buffer cells
                if (fluid.solid_mask[idx(fluid, i, j)] !== 0) { // only solve for fluid cells
                    
                    const div = (fluid.u[idx(fluid, i+1, j)] - fluid.u[idx(fluid, i, j)]) + // raw difference instead of proper divergenece
                                (fluid.v[idx(fluid, i, j+1)] - fluid.v[idx(fluid, i, j)]);  // which would require division by h
                    
                    div * fluid.density / dt; // convert to pressure units

                    const sx_left   = fluid.solid_mask[idx(fluid, i-1, j)]; // check if
                    const sx_right  = fluid.solid_mask[idx(fluid, i+1, j)]; // neighbouring
                    const sy_down   = fluid.solid_mask[idx(fluid, i, j-1)]; // cells are
                    const sy_up     = fluid.solid_mask[idx(fluid, i, j+1)]; // solid or fluid
                    
                    const s_coeff = sx_left + sx_right + sy_down + sy_up; // fluid neighbours (out of 4)
                    if (s_coeff === 0) continue; // skip if all neighbors are solid

                    let p = (-div / s_coeff) * over_relaxation; // average pressure contrib from neighbour cells + SOR
                    fluid.pressure[idx(fluid, i, j)] += p * fluid.density * fluid.h / dt; // convert to pressure units & accumulate over iters

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
export function applyBoundaryConditions(fluid: Fluid) {
    for (let i = 0; i < fluid.nx; i++) {
        fluid.u[idx(fluid, i, 0)] = fluid.u[idx(fluid, i, 1)]; // bottom b'ary
        fluid.u[idx(fluid, i, fluid.ny - 1)] = fluid.u[idx(fluid, i, fluid.ny - 2)]; // top b'ary
    }
    for (let j = 0; j < fluid.ny; j++) {
        fluid.v[idx(fluid, 0, j)] = fluid.v[idx(fluid, 1, j)]; // left b'ary
        fluid.v[idx(fluid, fluid.nx - 1, j)] = fluid.v[idx(fluid, fluid.nx - 2, j)]; // right b'ary
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
            dx = 0.5 * fluid.h; // center
            dy = 0.5 * fluid.h; // of cell
            p = fluid.smoke_density; break;
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

export function averageVelocity(fluid: Fluid, i: number, j: number, direction: String): number {
    if (direction === "u") {
        return 0.25 * (fluid.u[idx(fluid, i, j)] + fluid.u[idx(fluid, i+1, j)] +     // value of u where normally
                     fluid.u[idx(fluid, i, j-1)] + fluid.u[idx(fluid, i+1, j-1)])    // v is defined (horizontal face)
    } else {
        return 0.25 * (fluid.v[idx(fluid, i-1, j+1)] + fluid.v[idx(fluid, i, j+1)] + // value of v where normally
                         fluid.v[idx(fluid, i-1, j)] + fluid.v[idx(fluid, i, j)])    // u is defined (vertical face)
    }
}