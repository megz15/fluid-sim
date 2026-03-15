<script lang="ts">
    import * as lbm from "$lib/fluid_lbm";
    import { onMount } from "svelte";

    let canvas: HTMLCanvasElement;
    let animationFrame: number;
    let sim: lbm.FluidLBM;

    const nx = 200;
    const ny = 80;
    const tau = 0.53;
    const inletVelocity = 0.1;
    const stepsPerFrame = 15;

    onMount(() => {
        const canvasCtx = canvas.getContext("2d");
        if (!canvasCtx) {
            console.error("Canvas loading failure");
            return;
        }
        
        sim = lbm.initLBM(nx, ny, tau);
        const canvasPixels = canvasCtx.createImageData(nx, ny);
        
        // solid obstruction
        const cylinderRadius = 8;
        // break symmetry and force vortex shedding
        const cylinderCenter = { x: 40, y: 35 }; 
        
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                if ((i - cylinderCenter.x)**2 + (j - cylinderCenter.y)**2 <= cylinderRadius**2) {
                    sim.solid_mask[lbm.idx(sim, i, j)] = 0; // 0 = solid
                }
            }
        }

        function render() {
            // physics
            for (let s = 0; s < stepsPerFrame; s++) {
                lbm.step(sim, inletVelocity);
            }

            // render
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const grid_idx = lbm.idx(sim, i, j);
                    const pixel_idx = (j * nx + i) * 4;
                    
                    if (sim.solid_mask[grid_idx] === 0) {
                        // solid black cylinder
                        canvasPixels.data[pixel_idx] = 0;
                        canvasPixels.data[pixel_idx + 1] = 0;
                        canvasPixels.data[pixel_idx + 2] = 0;
                        canvasPixels.data[pixel_idx + 3] = 255;
                    } else {
                        // calc speed magnitude
                        const u = sim.u[grid_idx];
                        const v = sim.v[grid_idx];
                        const speed = Math.sqrt(u**2 + v**2);
                        
                        // scale speed to a 0-255 grayscale value
                        // divide by inletVelocity to normalize color scale, then boost contrast
                        const c = Math.max(0, Math.min(255, Math.floor((speed / inletVelocity) * 200)));
                        
                        canvasPixels.data[pixel_idx] = c;
                        canvasPixels.data[pixel_idx + 1] = c;
                        canvasPixels.data[pixel_idx + 2] = c;
                        canvasPixels.data[pixel_idx + 3] = 255;
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