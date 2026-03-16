<script lang="ts">
    import * as lbm from "$lib/fluid_lbm";
    import { onMount } from "svelte";

    let canvas: HTMLCanvasElement;
    let animation_frame: number;
    let fluid_domain: lbm.FluidLBM;

    const nx = 200;
    const ny = 80;
    const tau = 0.53;
    const inlet_velocity = 0.2;
    const steps_per_frame = 15;

    // obstruction
    let shape: String = "cylinder";
    const onstacle_size = 10;
    const cx = 40; const cy = 35;

    // RE = (inlet_velocity * characteristic_length) / nu
    $: nu = (tau - 0.5) / 3; // kinematic viscosity from LBM relaxation time
    $: re = (inlet_velocity * onstacle_size * 2) / nu;

    function resetSim() {
        if (!canvas) return;

        // todo: set density, gravity
        fluid_domain = lbm.initLBM(nx, ny, tau);

        // define obstacle
        for (let i = 1; i < nx - 1; i++) {
            for (let j = 1; j < ny - 1; j++) {
                if (shape === "cylinder") {
                    if ((i - cx)**2 + (j - cy)**2 <= onstacle_size**2) {
                        fluid_domain.solid_mask[lbm.idx(fluid_domain, i, j)] = 0;
                        // fluid_domain.u[lbm.idx(fluid_domain, i, j)] = 0; // Enforce no-slip
                    }
                } else if (shape === "square") {
                    if (Math.abs(i - cx) <= onstacle_size && Math.abs(j - cy) <= onstacle_size) {
                        fluid_domain.solid_mask[lbm.idx(fluid_domain, i, j)] = 0;
                        // fluid_domain.u[lbm.idx(fluid_domain, i, j)] = 0; // Enforce no-slip
                    }
                }
            }
        }
    }

    onMount(() => {
        resetSim();
        const canvas_ctx = canvas.getContext("2d");
        if (!canvas_ctx) return;
        const canvas_pixels = canvas_ctx.createImageData(nx, ny);

        function render() {
            // physics
            for (let s = 0; s < steps_per_frame; s++) {
                lbm.step(fluid_domain, inlet_velocity);
            }

            // render
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const grid_idx = lbm.idx(fluid_domain, i, j);
                    const pixel_idx = (j * nx + i) * 4;
                    
                    if (fluid_domain.solid_mask[grid_idx] === 0) {
                        // solid black cylinder
                        canvas_pixels.data[pixel_idx] = 0;
                        canvas_pixels.data[pixel_idx + 1] = 0;
                        canvas_pixels.data[pixel_idx + 2] = 0;
                        canvas_pixels.data[pixel_idx + 3] = 255;
                    } else {
                        // calc speed magnitude
                        const u = fluid_domain.u[grid_idx];
                        const v = fluid_domain.v[grid_idx];
                        const speed = Math.sqrt(u**2 + v**2);
                        const normalized_speed = Math.max(0, Math.min(1, speed / inlet_velocity));

                        canvas_pixels.data[pixel_idx] = Math.floor(normalized_speed * 100);
                        canvas_pixels.data[pixel_idx + 1] = Math.floor(normalized_speed * 255);
                        canvas_pixels.data[pixel_idx + 2] = Math.floor(200 + normalized_speed * 55);
                        canvas_pixels.data[pixel_idx + 3] = 255;
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
    <h2>Lattice Boltzmann Method (D2Q9)</h2>
    
    <div>
        <div>
            <label for="shape_lbm">Obstacle: </label>
            <select id="shape_lbm" bind:value={shape} on:change={resetSim}>
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