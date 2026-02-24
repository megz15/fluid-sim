export class FluidGrid {
    nx: number;
    ny: number;
    h: number;

    u: Float32Array;
    v: Float32Array;
    p: Float32Array;
    s: Float32Array;

    constructor(nx: number, ny: number, h: number) {
        this.nx = nx;
        this.ny = ny;
        this.h = h;

        const N = nx * ny;

        this.u = new Float32Array(N);
        this.v = new Float32Array(N);
        this.p = new Float32Array(N);
        this.s = new Float32Array(N);
    }

    idx(i: number, j: number): number {
        return i * this.ny + j;
    }
}