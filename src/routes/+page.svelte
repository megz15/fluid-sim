<script lang="ts">
    import * as fdm from "$lib/fluid_fdm";
    import * as lbm from "$lib/fluid_lbm";
    import { onMount } from "svelte";

    let canvas: HTMLCanvasElement;
    let animation_frame: number;

    let current_frame = 0;
    let v_history: number[] = [];
    let crossing_frames: number[] = [];
    let strouhal = 0;
    let karman_a = 0;
    let karman_b = 0;
    let karman_ratio = 0;

    // global params
    let method_used: "lbm" | "fdm" = "fdm";
    let wall_BC: "free-slip" | "no-slip" = "no-slip";
    const nx = 200;
    const ny = 80; // 2.5:1 tunnel
    // const gravity = 9.81; // m/s^2

    // colours
    let obstruction_colour = { r: 0, g: 0, b: 0 }; // black
    let flow_colour = {
        base: { r: 0, g: 0, b: 200 },
        variable: { r: 100, g: 255, b: 55 },
        a: 255
    };
    
    // obstruction
    let shape: String = "circle";
    let obstacle_size = 10;
    let cx = 40; let cy = 35;

    // eulerian FDM params
    let fdm_domain: fdm.Fluid;
    let fdm_nu = 0.001; // kinematic viscosity
    let fdm_h = 0.01; // grid spacing
    let fdm_dt = 1/60; // 60 FPS
    let fdm_inlet_velocity = 2; // m/s
    let fdm_diffuse_iters = 20; // number of diffusion iterations
    let fdm_pressure_iters = 40; // pressure solver iters
    let overrelaxation = 1.9; // SOR, 1.0 = Gauss-Seidel, <2.0 for convergence speedup

    // LBM D2Q9 params
    let lbm_domain: lbm.FluidLBM;
    let lbm_tau = 0.53; // relaxation time
    let lbm_inlet_velocity = 0.2; // lattice units
    let lbm_steps_per_frame = 15; // speed vs accuracy tradeoff
    let lbm_current_step = 0; // track total steps for inlet ramping
    let lbm_shock_ramp_steps = 1000; // steps to ramp up inl700et velocity to avoid "shock"

    // Profiling & Stability Metrics
    let compute_time_ms = 0; // milliseconds per frame
    let mlups = 0; // million lattice updates per second
    let max_speed = 0; // max velocity magnitude in the domain
    $: cfl = (max_speed * fdm_dt) / fdm_h; // Courant–Friedrichs–Lewy condition, <1
    $: mach = max_speed * Math.sqrt(3); // cs = 1/sqrt(3) in LBM

    // derived values
    // RE = (inlet_velocity * characteristic_length) / nu
    $: lbm_nu = (lbm_tau - 0.5) / 3;
    $: re_fdm = (fdm_inlet_velocity * obstacle_size * 2) * fdm_h / fdm_nu;
    $: re_lbm = (lbm_inlet_velocity * obstacle_size * 2) / lbm_nu;

    // reactivity
    $: {
        method_used, shape, obstacle_size, cx, cy, fdm_h, lbm_tau, lbm_nu, fdm_diffuse_iters, fdm_pressure_iters, lbm_shock_ramp_steps;
        if (canvas) resetSim();
    }

    function updateWakeMetrics() {
        current_frame++;
        
        const is_fdm = method_used === "fdm";
        const domain = is_fdm ? fdm_domain : lbm_domain as any;
        const lib = is_fdm ? fdm : lbm;
        
        // 1. Virtual Probe for Strouhal Number
        const probe_x = Math.floor(cx + obstacle_size * 4);
        const probe_y = Math.floor(cy + obstacle_size * 0.5);
        
        if (probe_x < nx && probe_y < ny) {
            const idx = lib.idx(domain, probe_x, probe_y);
            const current_v = domain.v[idx];
            v_history.push(current_v);
            
            // Detect zero-crossing (negative to positive v-velocity)
            if (v_history.length > 2) {
                const prev_v = v_history[v_history.length - 2];
                if (prev_v <= 0 && current_v > 0) {
                    crossing_frames.push(current_frame);
                    
                    if (crossing_frames.length > 2) {
                        let total_frames = 0;
                        for (let i = 1; i < crossing_frames.length; i++) {
                            total_frames += (crossing_frames[i] - crossing_frames[i-1]);
                        }
                        const avg_frames_per_cycle = total_frames / (crossing_frames.length - 1);
                        
                        const time_per_frame = is_fdm ? fdm_dt : lbm_steps_per_frame; 
                        const period = avg_frames_per_cycle * time_per_frame;
                        const frequency = 1.0 / period;
                        
                        const U = is_fdm ? fdm_inlet_velocity : lbm_inlet_velocity;
                        const D = obstacle_size * 2 * (is_fdm ? fdm_h : 1);
                        
                        strouhal = (frequency * D) / U;
                        
                        if (crossing_frames.length > 5) crossing_frames.shift();
                    }
                }
            }
            if (v_history.length > 100) v_history.shift();
        }

        // 2. Spatial Scanning for 'a' and 'b' parameters
        if (current_frame % 30 === 0) { 
            let max_vorticity = -Infinity;
            let min_vorticity = Infinity;
            let y_max_vort = cy;
            let y_min_vort = cy;
            
            const scan_x = probe_x;
            
            const scan_min = Math.max(1, cy - obstacle_size * 3);
            const scan_max = Math.min(ny - 2, cy + obstacle_size * 3);
            
            for (let j = scan_min; j <= scan_max; j++) {
                const v_right = domain.v[lib.idx(domain, scan_x + 1, j)];
                const v_left = domain.v[lib.idx(domain, scan_x - 1, j)];
                const u_up = domain.u[lib.idx(domain, scan_x, j + 1)];
                const u_down = domain.u[lib.idx(domain, scan_x, j - 1)];
                
                const dx = is_fdm ? (2 * fdm_h) : 2;
                const dy = is_fdm ? (2 * fdm_h) : 2;
                
                const vorticity = (v_right - v_left) / dx - (u_up - u_down) / dy;
                
                if (vorticity > max_vorticity) {
                    max_vorticity = vorticity;
                    y_max_vort = j;
                }
                if (vorticity < min_vorticity) {
                    min_vorticity = vorticity;
                    y_min_vort = j;
                }
            }
            
            karman_b = Math.abs(y_max_vort - y_min_vort) * (is_fdm ? fdm_h : 1);
            
            if (strouhal > 0) {
                const U = is_fdm ? fdm_inlet_velocity : lbm_inlet_velocity;
                const U_convection = 0.85 * U; 
                const frequency = (strouhal * U) / (obstacle_size * 2 * (is_fdm ? fdm_h : 1));
                karman_a = U_convection / frequency;
                
                karman_ratio = karman_b / karman_a;
            }
        }
    }

    function resetSim() {
        if (!canvas) return;
        lbm_current_step = 0;
        
        // Reset metrics
        current_frame = 0;
        v_history = [];
        crossing_frames = [];
        strouhal = 0;
        karman_a = 0;
        karman_b = 0;
        karman_ratio = 0;

        if (method_used === "fdm") {
            fdm_domain = fdm.initFluid(nx, ny, fdm_h);
        } else {
            lbm_domain = lbm.initLBM(nx, ny, lbm_tau);
        }

        const current_domain = method_used === "fdm" ? fdm_domain : lbm_domain;
        const lib = method_used === "fdm" ? fdm : lbm;

        for (let i = 1; i < nx - 1; i++) {
            for (let j = 1; j < ny - 1; j++) {
                const dx = i - cx;
                const dy = j - cy;
                let is_solid = false;

                if (shape === "circle") {
                    if (dx ** 2 + dy ** 2 <= obstacle_size ** 2) is_solid = true;

                } else if (shape === "square") {
                    if (Math.abs(dx) <= obstacle_size && Math.abs(dy) <= obstacle_size) is_solid = true;

                } else if (shape === "diamond") {
                    if (Math.abs(dx) + Math.abs(dy) <= obstacle_size * 1.2) is_solid = true;

                } else if (shape === "teardrop") {
                    if (dx >= 0) {
                        if (dx ** 2 + dy ** 2 <= (obstacle_size/1.5) ** 2) is_solid = true;
                    } else {
                        const tail_length = obstacle_size * 3;
                        if (dx >= -tail_length && Math.abs(dy) <= (obstacle_size/1.5) * (1 + dx / tail_length)) is_solid = true;
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
                    
                    if (method_used === "fdm") {
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

        const t0 = performance.now(); // Start timer

        if (method_used === "fdm") {
            for (let j = 0; j < ny; j++) {
                fdm_domain.u[fdm.idx(fdm_domain, 1, j)] = fdm_inlet_velocity;
                fdm_domain.v[fdm.idx(fdm_domain, 1, j)] = 0;
            }
            fdm.step(fdm_domain, fdm_dt, fdm_pressure_iters, overrelaxation, fdm_nu, fdm_diffuse_iters, wall_BC);
        } else {
            for (let s = 0; s < lbm_steps_per_frame; s++) {
                lbm_current_step++;
            
                // Ramp up over 1000 steps to avoid initial shock
                let ramp_factor = Math.min(1.0, lbm_current_step / lbm_shock_ramp_steps);

                const current_inlet_vel = lbm_inlet_velocity * ramp_factor;
                lbm.step(lbm_domain, current_inlet_vel, wall_BC);
            }
        }

        updateWakeMetrics();

        const t1 = performance.now(); // End timer
        compute_time_ms = t1 - t0;

        // Calculate MLUPS for LBM
        if (method_used === "lbm") {
            const total_updates = nx * ny * lbm_steps_per_frame;
            mlups = (total_updates / (compute_time_ms / 1000)) / 1e6;
        }

        // Reset max speed for this frame
        max_speed = 0;

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
        
        if (method_used === "fdm") {
            drawPixels(canvas_pixels, fdm_domain.solid_mask, fdm_domain.u, fdm_domain.v, fdm_inlet_velocity);
        } else {
            drawPixels(canvas_pixels, lbm_domain.solid_mask, lbm_domain.u, lbm_domain.v, lbm_inlet_velocity);
        }
        
        canvas_ctx.putImageData(canvas_pixels, 0, 0);
        animation_frame = requestAnimationFrame(render);
    }

    function drawPixels(canvas_pixels: ImageData, solid_mask: Float32Array, u: Float32Array, v: Float32Array, inlet_velocity: number) {
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                const grid_idx = method_used === "fdm"
                    ? (i + 1) * fdm_domain.ny + (j + 1)
                    : i * ny + j;
                const pixel_idx = (j * nx + i) * 4; // R, G, B, A channels
                
                if (solid_mask[grid_idx] === 0) {
                    canvas_pixels.data[pixel_idx] = obstruction_colour.r;
                    canvas_pixels.data[pixel_idx + 1] = obstruction_colour.g;
                    canvas_pixels.data[pixel_idx + 2] = obstruction_colour.b;
                    canvas_pixels.data[pixel_idx + 3] = 255; // black opaque
                } else {
                    // const smoke = fluid_domain.fluid_density[fluid_idx];
                    // const c = Math.max(0, Math.min(255, Math.floor(255 * smoke)));

                    // calculate speed
                    const speed = Math.sqrt(u[grid_idx]**2 + v[grid_idx]**2);
                    const normalized_speed = Math.max(0, Math.min(1, speed / inlet_velocity));

                    // Track max speed for CFL/Mach calculations
                    if (speed > max_speed) max_speed = speed;

                    canvas_pixels.data[pixel_idx] = Math.floor(flow_colour.base.r + normalized_speed * flow_colour.variable.r);
                    canvas_pixels.data[pixel_idx + 1] = Math.floor(flow_colour.base.g + normalized_speed * flow_colour.variable.g);
                    canvas_pixels.data[pixel_idx + 2] = Math.floor(flow_colour.base.b + normalized_speed * flow_colour.variable.b);
                    canvas_pixels.data[pixel_idx + 3] = flow_colour.a; // opaque
                }
            }
        }
    }

    onMount(() => {
        resetSim();
        render();
        return () => cancelAnimationFrame(animation_frame);
    });
</script>

<main class="p-6 max-w-7xl mx-auto min-h-screen font-sans text-slate-800">
    <header class="mb-6 border-b border-slate-200 pb-4 flex justify-between items-end">
        <div>
            <h1 class="text-3xl font-extrabold text-slate-900 tracking-tight">Computational Fluid Dynamics</h1>
            <p class="text-slate-500 text-sm">Real-time comparison of Eulerian FDM and Lattice Boltzmann Method</p>
        </div>
    </header>
    
    <div class="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-8">
        <section class="bg-white p-5 rounded-xl shadow-sm border border-slate-200">
            <h3 class="text-sm font-bold uppercase tracking-widest text-blue-600 mb-4 flex items-center gap-2">
                <span class="w-2 h-2 rounded-full bg-blue-600"></span> Global Config
            </h3>
            <div class="space-y-3">
                <label class="flex justify-between items-center text-sm">
                    <span class="font-medium">Solver Method</span>
                    <select bind:value={method_used} class="border rounded-md px-2 py-1 bg-slate-50 text-sm focus:ring-2 focus:ring-blue-500 outline-none">
                        <option value="fdm">Eulerian FDM</option>
                        <option value="lbm">LBM (D2Q9)</option>
                    </select>
                </label>
                <label class="flex justify-between items-center text-sm">
                    <span class="font-medium">Cylinder Geometry</span>
                    <select bind:value={shape} class="border rounded-md px-2 py-1 bg-slate-50 text-sm outline-none">
                        <option value="circle">Circular</option>
                        <option value="square">Square</option>
                        <option value="wedge">Triangular</option>
                        <option value="diamond">Diamond</option>
                        <option value="teardrop">Teardrop</option>
                        <option value="star">Star</option>
                        <option value="plate">Plate</option>
                    </select>
                </label>
                <label class="flex justify-between items-center text-sm">
                    <span class="font-medium">Boundary Condition</span>
                    <select bind:value={wall_BC} class="border rounded-md px-2 py-1 bg-slate-50 text-sm outline-none">
                        <option value="free-slip">Free-Slip</option>
                        <option value="no-slip">No-Slip</option>
                    </select>
                </label>
                <div class="pt-2 border-t border-slate-100 mt-2">
                    <div class="grid grid-cols-3 gap-3">
                        <label class="block">
                            <span class="text-[10px] uppercase font-bold text-slate-400">Obstacle Size</span>
                            <input type="number" bind:value={obstacle_size} min="2" max="40" class="w-full border rounded px-2 py-1 text-sm">
                        </label>
                        <label class="block">
                            <span class="text-[10px] uppercase font-bold text-slate-400">Position X</span>
                            <input type="number" bind:value={cx} min="10" max={nx-10} class="w-full border rounded px-2 py-1 text-sm">
                        </label>
                        <label class="block">
                            <span class="text-[10px] uppercase font-bold text-slate-400">Position Y</span>
                            <input type="number" bind:value={cy} min="10" max={ny-10} class="w-full border rounded px-2 py-1 text-sm">
                        </label>
                    </div>
                </div>
            </div>
        </section>

        <section class="bg-white p-5 rounded-xl shadow-sm border border-slate-200">
            {#if method_used === 'fdm'}
                <h3 class="text-sm font-bold uppercase tracking-widest text-emerald-600 mb-4 flex items-center gap-2">
                    <span class="w-2 h-2 rounded-full bg-emerald-600"></span> Finite Difference (Re: {re_fdm.toFixed(0)})
                </h3>
                <div class="space-y-2 text-sm">
                    <label class="flex justify-between items-center">Inlet Velocity <input type="number" step="0.1" bind:value={fdm_inlet_velocity} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label>
                    <label class="flex justify-between items-center">Viscosity (&nu;) <input type="number" step="0.0001" bind:value={fdm_nu} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label>
                    <!-- <label class="flex justify-between items-center">Grid Space (h) <input type="number" step="0.01" bind:value={fdm_h} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label> -->
                    <label class="flex justify-between items-center">Timestep (dt) <input type="number" step="0.001" bind:value={fdm_dt} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label>
                    <label class="flex justify-between items-center">Diffusion Iters <input type="number" bind:value={fdm_diffuse_iters} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label>
                    <label class="flex justify-between items-center">Pressure Iters <input type="number" bind:value={fdm_pressure_iters} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label>
                    <label class="flex justify-between items-center">SOR Overrelaxation <input type="number" step="0.01" min="0" max="5" bind:value={overrelaxation} class="w-20 border rounded px-2 py-0.5 text-right font-mono"></label>
                </div>
            {:else}
                <h3 class="text-sm font-bold uppercase tracking-widest text-purple-600 mb-4 flex items-center gap-2">
                    <span class="w-2 h-2 rounded-full bg-purple-600"></span> Lattice Boltzmann (Re: {re_lbm.toFixed(0)})
                </h3>
                <div class="space-y-2 text-sm">
                    <label class="flex justify-between items-center">Lattice Velocity <input type="number" step="0.01" bind:value={lbm_inlet_velocity} class="w-20 border rounded px-2 py-1 text-right font-mono"></label>
                    <label class="flex justify-between items-center">Relaxation (&tau;) <input type="number" step="0.01" min="0" bind:value={lbm_tau} class="w-20 border rounded px-2 py-1 text-right font-mono"></label>
                    <label class="flex justify-between items-center">Steps/Frame <input type="number" bind:value={lbm_steps_per_frame} class="w-20 border rounded px-2 py-1 text-right font-mono"></label>
                    <label class="flex justify-between items-center">Shock Ramp <input type="number" step="100" bind:value={lbm_shock_ramp_steps} class="w-20 border rounded px-2 py-1 text-right font-mono"></label>
                </div>
            {/if}
        </section>

        <section class="bg-white p-5 rounded-xl shadow-sm border border-slate-200">
            <h3 class="text-sm font-bold uppercase tracking-widest text-amber-600 mb-4 flex items-center gap-2">
                <span class="w-2 h-2 rounded-full bg-amber-600"></span> Colours
            </h3>
            <div class="space-y-4">
                <div>
                    <span class="text-[10px] uppercase font-extrabold text-slate-400 block mb-1">Obstruction RGB</span>
                    <div class="flex gap-1">
                        <input type="number" bind:value={obstruction_colour.r} min="0" max="255" class="w-full border rounded px-1 text-center bg-red-50 text-xs py-1" title="Red">
                        <input type="number" bind:value={obstruction_colour.g} min="0" max="255" class="w-full border rounded px-1 text-center bg-green-50 text-xs py-1" title="Green">
                        <input type="number" bind:value={obstruction_colour.b} min="0" max="255" class="w-full border rounded px-1 text-center bg-blue-50 text-xs py-1" title="Blue">
                    </div>
                </div>
                <div>
                    <span class="text-[10px] uppercase font-extrabold text-slate-400 block mb-1">Flow Base RGB</span>
                    <div class="flex gap-1">
                        <input type="number" bind:value={flow_colour.base.r} min="0" max="255" class="w-full border rounded px-1 text-center bg-red-50 text-xs py-1">
                        <input type="number" bind:value={flow_colour.base.g} min="0" max="255" class="w-full border rounded px-1 text-center bg-green-50 text-xs py-1">
                        <input type="number" bind:value={flow_colour.base.b} min="0" max="255" class="w-full border rounded px-1 text-center bg-blue-50 text-xs py-1">
                    </div>
                </div>
                <div>
                    <span class="text-[10px] uppercase font-extrabold text-slate-400 block mb-1">Flow Variable & Opacity</span>
                    <div class="flex gap-1">
                        <input type="number" bind:value={flow_colour.variable.r} min="0" max="255" class="w-full border rounded px-1 text-center bg-red-50 text-xs py-1" title="Red Speed Scaling">
                        <input type="number" bind:value={flow_colour.variable.g} min="0" max="255" class="w-full border rounded px-1 text-center bg-green-50 text-xs py-1" title="Green Speed Scaling">
                        <input type="number" bind:value={flow_colour.variable.b} min="0" max="255" class="w-full border rounded px-1 text-center bg-blue-50 text-xs py-1" title="Blue Speed Scaling">
                        <input type="number" bind:value={flow_colour.a} min="0" max="255" class="w-full border rounded px-1 text-center bg-white text-xs py-1 font-bold" title="Alpha (Opacity)">
                    </div>
                </div>
            </div>
        </section>
    </div>

    <div class="relative bg-slate-900 p-2 rounded-2xl shadow-2xl overflow-hidden border-8 border-slate-800">
        <canvas bind:this={canvas}
            width={nx} height={ny}
            class="w-full h-auto rounded-lg cursor-crosshair"
            style="image-rendering: pixelated;"
        ></canvas>
        <div class="absolute top-4 right-4 bg-slate-900/80 backdrop-blur text-white px-3 py-1 rounded text-[10px] font-mono uppercase tracking-widest border border-white/10">
            Method used: {method_used.toUpperCase()}
        </div>
        <div class="absolute bottom-4 right-4 bg-slate-900/80 backdrop-blur text-white px-3 py-1 rounded text-[10px] font-mono uppercase tracking-widest border border-white/10">
            {nx} x {ny} Lattice
        </div>
    </div>

    <div class="mt-8 grid grid-cols-2 md:grid-cols-5 gap-4">
        <div class="bg-white p-4 rounded-xl shadow-sm border border-slate-200 group">
            <span class="block text-[10px] uppercase font-bold text-slate-400 mb-1">Compute Time</span>
            <span class="text-xl font-mono font-bold text-slate-700">{compute_time_ms.toFixed(2)}<span class="text-xs ml-1 font-normal text-slate-400">ms</span></span>
        </div>
        
        {#if method_used === 'lbm'}
            <div class="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
                <span class="block text-[10px] uppercase font-bold text-slate-400 mb-1">Throughput</span>
                <span class="text-xl font-mono font-bold text-slate-700">{mlups.toFixed(2)}<span class="text-xs ml-1 font-normal text-slate-400">MLUPS</span></span>
            </div>
            <div class="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
                <span class="block text-[10px] uppercase font-bold text-slate-400 mb-1">Max Mach</span>
                <span class="text-xl font-mono font-bold {mach > 0.3 ? 'text-red-500' : 'text-slate-700'}">{mach.toFixed(3)}</span>
            </div>
        {:else}
            <div class="bg-white p-4 rounded-xl shadow-sm border border-slate-200 col-span-2">
                <span class="block text-[10px] uppercase font-bold text-slate-400 mb-1">CFL Stability Number</span>
                <span class="text-xl font-mono font-bold {cfl > 5 ? 'text-red-500' : 'text-slate-700'}">{cfl.toFixed(3)}</span>
            </div>
        {/if}

        <div class="bg-white p-4 rounded-xl shadow-sm border border-slate-200">
            <span class="block text-[10px] uppercase font-bold text-slate-400 mb-1">Strouhal (St)</span>
            <span class="text-xl font-mono font-bold text-blue-600">{strouhal > 0 ? strouhal.toFixed(3) : '---'}</span>
        </div>

        <div class="bg-white p-3 rounded-xl shadow-sm border border-slate-200 flex flex-col justify-between">
            <span class="block text-[10px] uppercase font-bold text-slate-400 mb-1">Karman Spacing</span>
            <div class="flex items-end justify-between">
                <span class="text-xl font-mono font-bold {Math.abs(karman_ratio - 0.281) < 0.05 ? 'text-emerald-600' : 'text-slate-700'}">
                    {karman_ratio > 0 ? karman_ratio.toFixed(3) : '---'} <span class="text-xs font-normal text-slate-400 ml-1">b/a</span>
                </span>
                <div class="text-[10px] font-mono text-slate-500 text-right leading-tight border-l border-slate-200 pl-2">
                    <span class="text-slate-400">a:</span> {karman_a > 0 ? karman_a.toFixed(2) : '--'}<br>
                    <span class="text-slate-400">b:</span> {karman_b > 0 ? karman_b.toFixed(2) : '--'}
                </div>
            </div>
        </div>
    </div>
</main>