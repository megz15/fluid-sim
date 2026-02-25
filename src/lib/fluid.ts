// FDM + Chorin's Projection
// MAC staggered grid

// todo: refactor monolithic code into separate files

export class Fluid {

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

    constructor(density: number, nx: number, ny: number, h: number) {
        this.density = density; // assuming constant density throughout
        this.nx = nx + 2; // adding buffer cells
        this.ny = ny + 2; // for boundary handling
        this.h = h; // cell size

        this.N = nx * ny; // total cells

        this.u = new Float32Array(this.N); // horizontal velocity
        this.v = new Float32Array(this.N); // vertical velocity
        this.pressure = new Float32Array(this.N);
        this.solid_mask = new Float32Array(this.N); // 0 for solid, 1 for fluid
        this.smoke_density = new Float32Array(this.N);

        this.smoke_density.fill(1); // clear fluid
    }

    idx(i: number, j: number): number {
        return i * this.ny + j;
    }

    // Euler ODE approximation
    integrateGravity(g: number, dt: number) {
        for (let i = 1; i < this.nx - 1; i++) {
            for (let j = 1; j < this.ny - 1; j++) { // leave out buffer cells
                if (this.solid_mask[this.idx(i, j)] !== 0 && this.solid_mask[this.idx(i, j-1)]) { // skip if cell b'ary above and below is solid
                    this.v[this.idx(i, j)] += g * dt; // apply gravity to vertical velocity
                }
            }
        }
    }

    // Incompressibility divergence-free condition
    solvePressure(iters: number, dt: number, over_relaxation: number) {
        for (let iter = 0; iter < iters; iter++) {
            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) { // leave out buffer cells
                    if (this.solid_mask[this.idx(i, j)] !== 0) { // only solve for fluid cells
                        
                        const div = (this.u[this.idx(i+1, j)] - this.u[this.idx(i, j)]) +   // raw difference instead of proper divergenece
                                    (this.v[this.idx(i, j+1)] - this.v[this.idx(i, j)]);    // which would require division by h
                        
                        div * this.density / dt; // convert to pressure units

                        const sx_left   = this.solid_mask[this.idx(i-1, j)]; // check if
                        const sx_right  = this.solid_mask[this.idx(i+1, j)]; // neighbouring
                        const sy_down   = this.solid_mask[this.idx(i, j-1)]; // cells are
                        const sy_up     = this.solid_mask[this.idx(i, j+1)]; // solid or fluid
                        
                        const s_coeff = sx_left + sx_right + sy_down + sy_up;
                        if (s_coeff === 0) continue; // skip if all neighbors are solid

                        let p = (-div / s_coeff) * over_relaxation * this.density * this.h / dt;
                        this.pressure[this.idx(i, j)] += p;

                        if (sx_left !== 0) this.u[this.idx(i, j)] -= p;
                        if (sx_right !== 0) this.u[this.idx(i+1, j)] += p; // apply pressure gradient to horizontal velocity
                        if (sy_down !== 0) this.v[this.idx(i, j)] -= p;
                        if (sy_up !== 0) this.v[this.idx(i, j+1)] += p; // apply pressure gradient to vertical velocity
                    }
                }
            }
        }
    }

    // free-slip (neumann, (u_i+1 - u_i)/h = 0) b'ary conditions (zero shear stress at b'ary)
    applyBoundaryConditions() {
        for (let i = 0; i < this.nx; i++) {
            this.u[this.idx(i, 0)] = this.u[this.idx(i, 1)]; // bottom b'ary
            this.u[this.idx(i, this.ny - 1)] = this.u[this.idx(i, this.ny - 2)]; // top b'ary
        }
        for (let j = 0; j < this.ny; j++) {
            this.v[this.idx(0, j)] = this.v[this.idx(1, j)]; // left b'ary
            this.v[this.idx(this.nx - 1, j)] = this.v[this.idx(this.nx - 2, j)]; // right b'ary
        }
    }

    // parameter value at arbitrary position
    arbitraryPosnParameter(x: number, y: number, param: String) {
        x = Math.max(this.h, Math.min(x, (this.nx - 2) * this.h)); // clamp to
        y = Math.max(this.h, Math.min(y, (this.ny - 2) * this.h)); // valid grid range

        let dx = 0, dy = 0, p;
        switch (param) {
            case "u":
                dy = 0.5 * this.h; // u is staggered in vertical face
                p = this.u; break;
            case "v":
                dx = 0.5 * this.h; // v is staggered in horizontal face
                p = this.v; break;
            default:
                dx = 0.5 * this.h; // center
                dy = 0.5 * this.h; // of cell
                p = this.smoke_density; break;
        }

        const x0 = Math.floor((x - dx) / this.h); // integer index of grid node to the left of x
        const y0 = Math.floor((y - dy) / this.h); // integer index of grid node below y
        
        const frac_x_from_left = (x - dx) / this.h - x0; // fractional distance of x from left grid node
        const frac_y_from_left = (y - dy) / this.h - y0; // fractional distance of y from bottom grid node
        
        const frac_x_from_right = 1 - frac_x_from_left; // fractional distance of x from right grid node
        const frac_y_from_right = 1 - frac_y_from_left; // fractional distance of y from top grid node
        
        const x1 = x0 + 1; // integer index of grid node to the right of x
        const y1 = y0 + 1; // integer index of grid node above y

        // bilinear interpolation
        return p[this.idx(x0, y0)] * frac_x_from_right * frac_y_from_right + // bottom left corner weighted by top right quad area
               p[this.idx(x1, y0)] * frac_x_from_left * frac_y_from_right +  // bottom right corner weighted by top left quad area
               p[this.idx(x0, y1)] * frac_x_from_right * frac_y_from_left +  // top left corner weighted by bottom right quad area
               p[this.idx(x1, y1)] * frac_x_from_left * frac_y_from_left;    // top right corner weighted by bottom left quad area
    }
}