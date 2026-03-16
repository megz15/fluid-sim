<!-- flow around infinitely long cylinder, with characteristic von karman vortex trails in wake region -->

<script lang="ts">
    import * as fluid from "$lib/fluid_eul";
    import { onMount } from "svelte";

    let canvas: HTMLCanvasElement;
    let animationFrame: number;
    let fluidDomain: fluid.Fluid;

    const density = 1000; // kg/m^3, water
    const gravity = 0; // m/s^2, no gravity for now
    const nu = 0.0015; // kinematic viscosity
    const nx = 200;
    const ny = 80; // 2.5:1 tunnel
    const h = 0.01;
    const dt = 1/60; // 60 FPS
    const inlet_velocity = 1.5; // m/s
    const diffuse_iters = 20; // number of diffusion iterations
    const pressure_iters = 40; // pressure solver iters
    const overrelaxation = 1.9;

    onMount(() => {
        const canvasCtx = canvas.getContext("2d");
        if (!canvasCtx) {
            console.error("Canvas loading failure");
            return;
        }
        
        fluidDomain = fluid.initFluid(density, nx, ny, h);
        const canvasPixels = canvasCtx.createImageData(nx, ny);
        
        // define cylinder
        const cylinderRadius = 8;
        const cylinderCenter = { x: 40, y: 35 };
        for (let i = 1; i < nx - 1; i++) {
            for (let j = 1; j < ny - 1; j++) {
                // (x-a)^2 + (y-b)^2 <= r^2
                if ((i - cylinderCenter.x)**2 + (j - cylinderCenter.y)**2 <= cylinderRadius**2) {
                    fluidDomain.solid_mask[fluid.idx(fluidDomain, i, j)] = 0;
                }
            }
        }

        // fluidDomain.u[fluid.idx(fluidDomain, 100, 20)] = 20; 

        function render() {
            for (let j = 0; j < ny; j++) {
                fluidDomain.u[fluid.idx(fluidDomain, 1, j)] = 1.5;
                // fluidDomain.fluid_density[fluid.idx(fluidDomain, 1, j)] = 0;
                if (j % 10 < 5) {
                    fluidDomain.fluid_density[fluid.idx(fluidDomain, 1, j)] = 0.0; // Black
                } else {
                    fluidDomain.fluid_density[fluid.idx(fluidDomain, 1, j)] = 1.0; // Clear
                }
            }

            fluid.step(fluidDomain, gravity, dt, pressure_iters, overrelaxation, nu, diffuse_iters);

            // render on canvas
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    const fluid_idx = fluid.idx(fluidDomain, i + 1, j + 1);
                    const pixel_idx = (j * nx + i) * 4; // R, G, B, A channels
                    
                    if (fluidDomain.solid_mask[fluid_idx] === 0) {
                        canvasPixels.data[pixel_idx] = 0;
                        canvasPixels.data[pixel_idx + 1] = 0;
                        canvasPixels.data[pixel_idx + 2] = 0;
                        canvasPixels.data[pixel_idx + 3] = 255; // black opaque
                    } else {
                        const smoke = fluidDomain.fluid_density[fluid_idx];
                        const c = Math.max(0, Math.min(255, Math.floor(255 * smoke)));
                        canvasPixels.data[pixel_idx] = c+40;
                        canvasPixels.data[pixel_idx + 1] = c+60;
                        canvasPixels.data[pixel_idx + 2] = c+100;
                        canvasPixels.data[pixel_idx + 3] = 255; // grey opaque
                    }
                }
            }
            canvasCtx!.putImageData(canvasPixels, 0, 0);
            animationFrame = requestAnimationFrame(render);
        }

        render();
        return () => {
            cancelAnimationFrame(animationFrame);
        };
    });
</script>

<main>
    <canvas bind:this={canvas}
        width={nx} height={ny}
        class="w-200 h-80 bg-white border-2"
        style="image-rendering: pixelated;"
    ></canvas>
</main>