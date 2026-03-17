<script lang="ts">
    import * as eul from "$lib/fluid_eul";
    import * as lbm from "$lib/fluid_lbm";
    import { onMount, onDestroy } from "svelte";

    let canvas: HTMLCanvasElement;
    let animation_frame: number;

    // global params
    let method_used: "lbm" | "eul" = "eul";
    const nx = 200;
    const ny = 80; // 2.5:1 tunnel
    const gravity = 0; // m/s^2, no gravity for now
    const overrelaxation = 1.9;
    
    // obstruction
    let shape: String = "cylinder";
    let obstacle_size = 10;
    let cx = 40; let cy = 35;

    // eulerian FDM params
    let eul_domain: eul.Fluid;
    let eul_density = 1000; // kg/m^3, water
    let eul_nu = 0.001; // kinematic viscosity
    let eul_h = 0.01;
    let eul_dt = 1/60; // 60 FPS
    let eul_inlet_velocity = 2; // m/s
    let eul_diffuse_iters = 20; // number of diffusion iterations
    let eul_pressure_iters = 40; // pressure solver iters

    // LBM D2Q9 params
    let lbm_domain: lbm.FluidLBM;
    let lbm_tau = 0.53; // relaxation time
    let lbm_inlet_velocity = 0.2; // lattice units
    let lbm_steps_per_frame = 15;

    // derived values
    // RE = (inlet_velocity * characteristic_length) / nu
    $: lbm_nu = (lbm_tau - 0.5) / 3;
    $: re_eul = (eul_inlet_velocity * obstacle_size * 2) / eul_nu;
    $: re_lbm = (lbm_inlet_velocity * obstacle_size * 2) / lbm_nu;

    // reactivity
    $: {
        method_used, shape, obstacle_size, cx, cy, eul_nu, lbm_tau, lbm_nu;
        if (canvas) resetSim();
    }

    function resetSim() {
        if (!canvas) return;

        // Initialize the appropriate domain
        if (method_used === "eul") {
            eul_domain = eul.initFluid(eul_density, nx, ny, eul_h);
        } else {
            lbm_domain = lbm.initLBM(nx, ny, lbm_tau);
        }

        const current_domain = method_used === "eul" ? eul_domain : lbm_domain;
        const lib = method_used === "eul" ? eul : lbm;

        // Single loop for obstacle definition
        for (let i = 1; i < nx - 1; i++) {
            for (let j = 1; j < ny - 1; j++) {
                const dx = i - cx;
                const dy = j - cy;
                let is_solid = false;

                if (shape === "cylinder") {
                    if (dx ** 2 + dy ** 2 <= obstacle_size ** 2) is_solid = true;

                } else if (shape === "square") {
                    if (Math.abs(dx) <= obstacle_size && Math.abs(dy) <= obstacle_size) is_solid = true;

                } else if (shape === "diamond") {
                    if (Math.abs(dx) + Math.abs(dy) <= obstacle_size * 1.2) is_solid = true;

                } else if (shape === "teardrop") {
                    if (dx <= 0) {
                        if (dx ** 2 + dy ** 2 <= (obstacle_size/1.5) ** 2) is_solid = true;
                    } else {
                        const tail_length = obstacle_size * 3;
                        if (dx <= tail_length && Math.abs(dy) <= (obstacle_size/1.5) * (1 - dx / tail_length)) is_solid = true;
                    }
                } else if (shape === "star") {
                    const angle = Math.atan2(dy, dx);
                    const r = Math.sqrt(dx ** 2 + dy ** 2);
                    // 5-point star formula
                    const r_bound = obstacle_size * (0.6 + 0.4 * Math.cos(5 * angle));
                    if (r <= r_bound) is_solid = true;

                } else if (shape === "wedge") {
                    // Triangle pointing left into the flow
                    if (dx >= -obstacle_size && dx <= obstacle_size) {
                        if (Math.abs(dy) <= (dx + obstacle_size) * 0.5) is_solid = true;
                    }

                } else if (shape === "plate") {
                    // Thin plate at a 45-degree angle
                    const angle = Math.PI / 4;
                    const rx = dx * Math.cos(angle) - dy * Math.sin(angle);
                    const ry = dx * Math.sin(angle) + dy * Math.cos(angle);
                    if (Math.abs(rx) <= obstacle_size * 0.15 && Math.abs(ry) <= obstacle_size * 1.5) is_solid = true;
                }

                if (is_solid) {
                    const active_domain = current_domain as any;

                    const index = lib.idx(active_domain, i, j);
                    current_domain.solid_mask[index] = 0;
                    
                    // For Eulerian, we explicitly zero out the velocity inside solids
                    if (method_used === "eul") {
                        current_domain.u[index] = 0;
                        current_domain.v[index] = 0;
                    }
                }
            }
        }
    }

    function render() {
        if (!canvas) return;
        const canvas_ctx = canvas.getContext("2d");
        if (!canvas_ctx) return;
        const canvas_pixels = canvas_ctx!.createImageData(nx, ny);

        if (method_used === "eul") {
            for (let j = 0; j < ny; j++) {
                eul_domain.u[eul.idx(eul_domain, 1, j)] = eul_inlet_velocity;
                eul_domain.v[eul.idx(eul_domain, 1, j)] = 0;
            }
            
            eul.step(eul_domain, gravity, eul_dt, eul_pressure_iters, overrelaxation, eul_nu, eul_diffuse_iters);
            
            drawPixels(canvas_pixels, eul_domain.solid_mask, eul_domain.u, eul_domain.v, eul_inlet_velocity);
        } else {
            for (let s = 0; s < lbm_steps_per_frame; s++) {
                lbm.step(lbm_domain, lbm_inlet_velocity);
            }
            drawPixels(canvas_pixels, lbm_domain.solid_mask, lbm_domain.u, lbm_domain.v, lbm_inlet_velocity);
        }

        // for (let j = 0; j < ny; j++) {
            // fluid_domain.u[fluid.idx(fluid_domain, 1, j)] = inlet_velocity;
            // fluid_domain.fluid_density[fluid.idx(fluid_domain, 1, j)] = 0;
            // if (j % 10 < 5) {
            //     fluid_domain.fluid_density[fluid.idx(fluid_domain, 1, j)] = 0.0; // Black
            // } else {
            //     fluid_domain.fluid_density[fluid.idx(fluid_domain, 1, j)] = 1.0; // Clear
            // }
        // }

        // render velocity magnitude
        canvas_ctx.putImageData(canvas_pixels, 0, 0);
        animation_frame = requestAnimationFrame(render);
    }

    function drawPixels(canvas_pixels: ImageData, solid_mask: Float32Array, u: Float32Array, v: Float32Array, inlet_velocity: number) {
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                const grid_idx = i * (ny+ (method_used === "eul" ? 2 : 0) ) + j;
                const pixel_idx = (j * nx + i) * 4; // R, G, B, A channels
                
                if (solid_mask[grid_idx] === 0) {
                    canvas_pixels.data[pixel_idx] = 0;
                    canvas_pixels.data[pixel_idx + 1] = 0;
                    canvas_pixels.data[pixel_idx + 2] = 0;
                    canvas_pixels.data[pixel_idx + 3] = 255; // black opaque
                } else {
                    // const smoke = fluid_domain.fluid_density[fluid_idx];
                    // const c = Math.max(0, Math.min(255, Math.floor(255 * smoke)));

                    // calculate speed
                    const speed = Math.sqrt(u[grid_idx]**2 + v[grid_idx]**2);
                    const normalized_speed = Math.max(0, Math.min(1, speed / inlet_velocity));

                    canvas_pixels.data[pixel_idx] = Math.floor(normalized_speed * 100);
                    canvas_pixels.data[pixel_idx + 1] = Math.floor(normalized_speed * 255);
                    canvas_pixels.data[pixel_idx + 2] = Math.floor(200 + normalized_speed * 55);
                    canvas_pixels.data[pixel_idx + 3] = 255; // opaque
                }
            }
        }
    }

    onMount(() => {
        resetSim();
        render();
        return () => cancelAnimationFrame(animation_frame);
    });

    onDestroy(() => {
        if (animation_frame) cancelAnimationFrame(animation_frame);
    });
</script>

<main>
    <h2>Fluid Dynamics Comparison</h2>
    
    <div>
        <div>
            <h3>Global Settings</h3>
            <label>Method: 
                <select bind:value={method_used}>
                    <option value="eul">Eulerian FDM</option>
                    <option value="lbm">Lattice Boltzmann (D2Q9)</option>
                </select>
            </label><br>
            <label>Obstacle Shape: 
                <select bind:value={shape}>
                    <option value="cylinder">Cylinder</option>
                    <option value="square">Square</option>
                    <option value="teardrop">Teardrop</option>
                    <option value="diamond">Diamond</option>
                    <option value="star">Star</option>
                    <option value="wedge">Wedge</option>
                    <option value="plate">Plate</option>
                </select>
            </label><br>
            <label>Obstacle Size: 
                <input type="number" bind:value={obstacle_size} min="2" max="30">
            </label>
        </div>

        {#if method_used === 'eul'}
            <div>
                <h3>Eulerian Params (Re: {re_eul.toFixed(1)})</h3>
                <label>Inlet Velocity: <input type="number" step="0.1" bind:value={eul_inlet_velocity}></label><br>
                <label>Kinematic Visc (&nu;): <input type="number" step="0.001" bind:value={eul_nu}></label><br>
                <label>Grid Spacing (h): <input type="number" step="0.01" bind:value={eul_h}></label><br>
                <label>Timestep (dt): <input type="number" step="0.001" bind:value={eul_dt}></label><br>
            </div>
        {:else}
            <div>
                <h3>LBM Params (Re: {re_lbm.toFixed(1)})</h3>
                <label>Inlet Velocity (lattice): <input type="number" step="0.01" bind:value={lbm_inlet_velocity}></label><br>
                <label>Relaxation Time (&tau;): <input type="number" step="0.01" min="0.51" bind:value={lbm_tau}></label><br>
                <label>Steps per Frame: <input type="number" step="1" bind:value={lbm_steps_per_frame}></label>
            </div>
        {/if}
    </div>

    <canvas bind:this={canvas}
        width={nx} height={ny}
        class="w-200 h-80 bg-white border-2"
        style="image-rendering: pixelated;"
    ></canvas>
</main>