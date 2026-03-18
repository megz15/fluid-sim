<script lang="ts">
	import './layout.css';
	import favicon from '$lib/assets/water-96.png';

	let { children } = $props();
</script>

<svelte:head>
	<link rel="icon" href={favicon} />
	<title>2D Fluid Simulation</title>
</svelte:head>
{@render children()}
