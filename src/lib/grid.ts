export class Fluid {

    // Chorin's Projection Method
    // MAC staggered grid

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
        this.density = density;
        this.nx = nx + 2; // adding buffer cells
        this.ny = ny + 2; // for boundary handling
        this.h = h;

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
                if (this.solid_mask[this.idx(i, j)] !== 0 && this.solid_mask[this.idx(i, j-1)]) { // skip if cell boundary above and below is solid
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

                        const sx_left = this.solid_mask[this.idx(i-1, j)];
                        const sx_right = this.solid_mask[this.idx(i+1, j)];
                        const sy_down = this.solid_mask[this.idx(i, j-1)];
                        const sy_up = this.solid_mask[this.idx(i, j+1)];
                        
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
}