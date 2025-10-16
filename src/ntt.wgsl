// Note: The shader is only valid  to be used with these values!
const n : i32 = 4;
const n_inverse : i32 = 5761;

struct wrapped_i32 {
  @size(16) elem: i32 // or `@align`...
}

const workgroup_size : i32 = 32;

@group(0) @binding(0) var<uniform> twiddles : array<wrapped_i32, n>;
@group(0) @binding(1) var<storage, read> input_a : array<i32, n>;
@group(0) @binding(2) var<storage, read> input_b: array<i32, n>;
@group(0) @binding(3) var<storage, read_write> output : array<i32, n>;

@compute @workgroup_size(workgroup_size)
fn add_twiddles(@builtin(global_invocation_id) id : vec3<u32>) {
  output[id.x] = input_a[id.x] + twiddles[id.x].elem;
}

@compute @workgroup_size(workgroup_size)
fn fwd(@builtin(global_invocation_id) id: vec3<u32>) {
  output[id.x] = input_a[id.x];
}
