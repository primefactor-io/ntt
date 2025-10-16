// See: https://gpuweb.github.io/gpuweb/wgsl/#address-space-layout-constraints
// See: https://eliemichel.github.io/LearnWebGPU/basic-3d-rendering/shader-uniforms/multiple-uniforms.html#memory-layout-constraints
struct wrapped_f32 {
  @size(16) elem: f32 // or `@align`...
}

const n : u32 = {n};

@group(0) @binding(0) var<uniform> constsBuffer : array<wrapped_f32, 64>;
@group(0) @binding(1) var<storage, read> inputBuffer : array<f32, 64>;
@group(0) @binding(2) var<storage, read_write> outputBuffer : array<f32, 64>;

fn f(x : f32) -> f32 {
  return (2.0 * x) + 1.0;
}

@compute @workgroup_size(32)
fn main(@builtin(global_invocation_id) id : vec3<u32>) {
  outputBuffer[id.x] = f(inputBuffer[id.x]) + constsBuffer[id.x].elem;
}
