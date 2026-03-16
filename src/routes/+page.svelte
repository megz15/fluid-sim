<!-- flow around infinitely long cylinder, with characteristic von karman vortex trails in wake region -->

<script lang="ts">
    import * as fluid from "$lib/fluid_eul";
    import { onMount } from "svelte";

    let canvas: HTMLCanvasElement;
    let animation_frame: number;
    let fluid_domain: fluid.Fluid;

    const density = 1000; // kg/m^3, water
    const gravity = 0; // m/s^2, no gravity for now
    const nu = 0.001; // kinematic viscosity
    const nx = 200;
    const ny = 80; // 2.5:1 tunnel
    const h = 0.01;
    const dt = 1/60; // 60 FPS
    const inlet_velocity = 1.5; // m/s
    const diffuse_iters = 20; // number of diffusion iterations
    const pressure_iters = 40; // pressure solver iters
    const overrelaxation = 1.9;

    // obstruction
    let shape: String = "cylinder";
    const onstacle_size = 10;
    const cx = 40; const cy = 35;

    // RE = (inlet_velocity * characteristic_length) / nu
    $: re = (inlet_velocity * onstacle_size * 2) / nu;

    function resetSim() {
        if (!canvas) return;
        const canvas_ctx = canvas.getContext("2d");
        if (!canvas_ctx) return;

        fluid_domain = fluid.initFluid(density, nx, ny, h);

        // todo: fix for init apparent "bulge" in sim

        // define obstacle
        for (let i = 1; i < nx - 1; i++) {
            for (let j = 1; j < ny - 1; j++) {
                if (shape === "cylinder") {
                    if ((i - cx)**2 + (j - cy)**2 <= onstacle_size**2) {
                        fluid_domain.solid_mask[fluid.idx(fluid_domain, i, j)] = 0;
                        fluid_domain.u[fluid.idx(fluid_domain, i, j)] = 0; // Enforce no-slip
                    }
                } else if (shape === "square") {
                    if (Math.abs(i - cx) <= onstacle_size && Math.abs(j - cy) <= onstacle_size) {
                        fluid_domain.solid_mask[fluid.idx(fluid_domain, i, j)] = 0;
                        fluid_domain.u[fluid.idx(fluid_domain, i, j)] = 0; // Enforce no-slip
                    }
                }
            }
        }
    }

    onMount(() => {
        resetSim();
        const canvas_ctx = canvas.getContext("2d");
        const canvas_pixels = canvas_ctx!.createImageData(nx, ny);

        // fluid_domain.u[fluid.idx(fluid_domain, 100, 20)] = 20; 

        function render() {
            for (let j = 0; j < ny; j++) {
                fluid_domain.u[fluid.idx(fluid_domain, 1, j)] = inlet_velocity;
                fluid_domain.v[fluid.idx(fluid_domain, 1, j)] = 0;
                // fluid_domain.fluid_density[fluid.idx(fluid_domain, 1, j)] = 0;
                // if (j % 10 < 5) {
                //     fluid_domain.fluid_density[fluid.idx(fluid_domain, 1, j)] = 0.0; // Black
                // } else {
                //     fluid_domain.fluid_density[fluid.idx(fluid_domain, 1, j)] = 1.0; // Clear
                // }
            }

            fluid.step(fluid_domain, gravity, dt, pressure_iters, overrelaxation, nu, diffuse_iters);

            // render velocity magnitude
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    const fluid_idx = fluid.idx(fluid_domain, i + 1, j + 1);
                    const pixel_idx = (j * nx + i) * 4; // R, G, B, A channels
                    
                    if (fluid_domain.solid_mask[fluid_idx] === 0) {
                        canvas_pixels.data[pixel_idx] = 0;
                        canvas_pixels.data[pixel_idx + 1] = 0;
                        canvas_pixels.data[pixel_idx + 2] = 0;
                        canvas_pixels.data[pixel_idx + 3] = 255; // black opaque
                    } else {
                        // const smoke = fluid_domain.fluid_density[fluid_idx];
                        // const c = Math.max(0, Math.min(255, Math.floor(255 * smoke)));

                        // calculate speed
                        const u = fluid_domain.u[fluid_idx];
                        const v = fluid_domain.v[fluid_idx];
                        const speed = Math.sqrt(u**2 + v**2);
                        const normalized_speed = Math.max(0, Math.min(1, speed / inlet_velocity));

                        canvas_pixels.data[pixel_idx] = Math.floor(normalized_speed * 100);
                        canvas_pixels.data[pixel_idx + 1] = Math.floor(normalized_speed * 255);
                        canvas_pixels.data[pixel_idx + 2] = Math.floor(200 + normalized_speed * 55);
                        canvas_pixels.data[pixel_idx + 3] = 255; // opaque
                    }
                }
            }
            canvas_ctx!.putImageData(canvas_pixels, 0, 0);
            animation_frame = requestAnimationFrame(render);
        }

        render();
        return () => cancelAnimationFrame(animation_frame);
    });
</script>

<main>
    <h2>Eulerian FDM</h2>
    
    <div>
        <div>
            <label for="shape">Obstacle: </label>
            <select id="shape" bind:value={shape} on:change={resetSim}>
                <option value="cylinder">Cylinder</option>
                <option value="square">Square</option>
            </select>
        </div>
        <div>
            <strong>Reynolds Number (Re):</strong> {re}
        </div>
    </div>

    <canvas bind:this={canvas}
        width={nx} height={ny}
        class="w-200 h-80 bg-white border-2"
        style="image-rendering: pixelated;"
    ></canvas>
</main>