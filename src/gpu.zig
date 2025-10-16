const std = @import("std");
const c = @cImport({
    @cInclude("webgpu/wgpu.h");
    @cInclude("webgpu/webgpu.h");
});

// const shader_code =
//     \\ // See: https://gpuweb.github.io/gpuweb/wgsl/#address-space-layout-constraints
//     \\ // See: https://eliemichel.github.io/LearnWebGPU/basic-3d-rendering/shader-uniforms/multiple-uniforms.html#memory-layout-constraints
//     \\ struct wrapped_f32 {
//     \\   @size(16) elem: f32 // or `@align`...
//     \\ }
//     \\
//     \\ const n : i32 = {n};
//     \\
//     \\ @group(0) @binding(0) var<uniform> constsBuffer : array<wrapped_f32, n>;
//     \\ @group(0) @binding(1) var<storage, read> inputBuffer : array<f32, n>;
//     \\ @group(0) @binding(2) var<storage, read_write> outputBuffer : array<f32, n>;
//     \\ fn f(x : f32) -> f32 {
//     \\   return (2.0 * x) + 1.0;
//     \\ }
//     \\
//     \\ @compute @workgroup_size(32)
//     \\ fn main(@builtin(global_invocation_id) id : vec3<u32>) {
//     \\   outputBuffer[id.x] = f(inputBuffer[id.x]) + constsBuffer[id.x].elem;
//     \\ }
// ;

pub fn Application(comptime buffer_size: u32) type {
    return struct {
        const Self = @This();

        comptime buffer_size: u32 = buffer_size,
        device: c.WGPUDevice,
        queue: c.WGPUQueue,
        bind_group_layout: c.WGPUBindGroupLayout,
        pipeline_layout: c.WGPUPipelineLayout,
        pipeline: c.WGPUComputePipeline,
        consts_buffer: c.WGPUBuffer,
        input_buffer: c.WGPUBuffer,
        output_buffer: c.WGPUBuffer,
        map_buffer: c.WGPUBuffer,
        bind_group: c.WGPUBindGroup,

        pub fn init() !Self {
            const device = try initDevice();
            const queue = c.wgpuDeviceGetQueue(device);

            const bind_group_layout = getBindGroupLayout(device);
            const compute_pipeline = getComputePipeline(device, bind_group_layout);

            const buffers = getBuffers(device, buffer_size);
            const bind_group = getBindGroup(device, bind_group_layout, buffer_size, buffers.consts, buffers.input, buffers.output);

            return Self{
                .buffer_size = buffer_size,
                .device = device,
                .queue = queue,
                .bind_group_layout = bind_group_layout,
                .pipeline_layout = compute_pipeline.layout,
                .pipeline = compute_pipeline.pipeline,
                .consts_buffer = buffers.consts,
                .input_buffer = buffers.input,
                .output_buffer = buffers.output,
                .map_buffer = buffers.map,
                .bind_group = bind_group,
            };
        }

        pub fn deinit(self: Self) void {
            c.wgpuBindGroupRelease(self.bind_group);
            c.wgpuBufferDestroy(self.consts_buffer);
            c.wgpuBufferRelease(self.consts_buffer);
            c.wgpuBufferDestroy(self.input_buffer);
            c.wgpuBufferRelease(self.input_buffer);
            c.wgpuBufferDestroy(self.output_buffer);
            c.wgpuBufferRelease(self.output_buffer);
            c.wgpuBufferDestroy(self.map_buffer);
            c.wgpuBufferRelease(self.map_buffer);
            c.wgpuComputePipelineRelease(self.pipeline);
            c.wgpuPipelineLayoutRelease(self.pipeline_layout);
            c.wgpuBindGroupLayoutRelease(self.bind_group_layout);
            c.wgpuQueueRelease(self.queue);
            c.wgpuDeviceRelease(self.device);
        }

        pub fn run(self: Self) void {
            const num_elements: usize = self.buffer_size / @sizeOf(f32);

            // Populate constants (uniform buffer).
            var consts = [_]f32{0} ** num_elements;
            for (0..num_elements) |i| {
                consts[i] = @as(f32, @floatFromInt(i));
            }

            c.wgpuQueueWriteBuffer(self.queue, self.consts_buffer, 0, &consts, consts.len * @sizeOf(f32));

            // Populate input.
            var input = [_]f32{0} ** num_elements;

            for (0..num_elements) |i| {
                input[i] = 0.1 * @as(f32, @floatFromInt(i));
            }

            c.wgpuQueueWriteBuffer(self.queue, self.input_buffer, 0, &input, input.len * @sizeOf(f32));

            // Initialize command encoder.
            const command_encoder_desc = c.WGPUCommandEncoderDescriptor{
                .label = "Command Encoder",
                .nextInChain = null,
            };
            const command_encoder = c.wgpuDeviceCreateCommandEncoder(self.device, &command_encoder_desc);
            defer c.wgpuCommandEncoderRelease(command_encoder);

            // Create compute pass.
            const compute_pass_desc = c.WGPUComputePassDescriptor{
                .label = "Compute Pass",
                .nextInChain = null,
                .timestampWrites = null,
            };
            const compute_pass = c.wgpuCommandEncoderBeginComputePass(command_encoder, &compute_pass_desc);
            defer c.wgpuComputePassEncoderRelease(compute_pass);

            // Use compute pass.
            c.wgpuComputePassEncoderSetPipeline(compute_pass, self.pipeline);
            c.wgpuComputePassEncoderSetBindGroup(compute_pass, 0, self.bind_group, 0, null);

            const invocation_count = self.buffer_size / @sizeOf(f32);
            const workgroup_size = 32;
            // ceil(invocation_count / workgroup_size)
            const workgroup_count = ((invocation_count + workgroup_size) - 1) / workgroup_size;
            c.wgpuComputePassEncoderDispatchWorkgroups(compute_pass, workgroup_count, 1, 2);

            // Finalize compute pass.
            c.wgpuComputePassEncoderEnd(compute_pass);

            // Before encoder finished.
            c.wgpuCommandEncoderCopyBufferToBuffer(command_encoder, self.output_buffer, 0, self.map_buffer, 0, self.buffer_size);

            // Encode and submit the GPU commands.
            const command_buffer_desc = c.WGPUCommandBufferDescriptor{
                .label = "Command Buffer",
                .nextInChain = null,
            };
            const commands = c.wgpuCommandEncoderFinish(command_encoder, &command_buffer_desc);
            defer c.wgpuCommandBufferRelease(commands);

            c.wgpuQueueSubmit(self.queue, 1, &commands);

            var buffer_user_data = BufferMapCallbackData{
                .buffer = self.map_buffer,
                .buffer_size = self.buffer_size,
            };

            c.wgpuBufferMapAsync(self.map_buffer, c.WGPUMapMode_Read, 0, self.buffer_size, bufferMapCallback, &buffer_user_data);

            while (buffer_user_data.output == null) {
                c.wgpuQueueSubmit(self.queue, 0, null); // Mimmicking the `processEvent` function as wgpu-native doesn't come with it yet.
                _ = c.wgpuDevicePoll(self.device, 0, null); // TODO: Do we need this? If not, we can remove the include at the top.
            }

            // const output = buffer_user_data.output orelse @panic("Output value not available");

            // std.debug.print("Output: {d}\n", .{output});
        }
    };
}

fn initDevice() !c.WGPUDevice {
    // Create WebGPU instance.
    const instance_desc = c.WGPUInstanceDescriptor{ .nextInChain = null };
    const instance = c.wgpuCreateInstance(&instance_desc) orelse {
        return error.CreateWebGPUInstance;
    };
    // The instance isn't needed once the adapter is selected.
    // See: https://eliemichel.github.io/LearnWebGPU/getting-started/adapter-and-device/the-adapter.html#destruction
    defer c.wgpuInstanceRelease(instance);

    // Create WebGPU adapter.
    const adapter_opts = c.WGPURequestAdapterOptions{ .nextInChain = null };
    var adapter_user_data = RequestAdapterCallbackData{};
    // When this function returns we know that our callback function was called (there's no need to await it).
    c.wgpuInstanceRequestAdapter(instance, &adapter_opts, requestAdapterCallback, &adapter_user_data);

    if (adapter_user_data.adapter == null) {
        return error.RequestWebGPUAdapter;
    }

    // The adapter isn't needed once the device is selected.
    // See: https://eliemichel.github.io/LearnWebGPU/getting-started/adapter-and-device/the-device.html#device-request
    const adapter = adapter_user_data.adapter.?;
    defer c.wgpuAdapterRelease(adapter);

    // Create WebGPU device.
    const device_desc = c.WGPUDeviceDescriptor{
        .label = "Compute Device",
        .nextInChain = null,
        .requiredFeatureCount = 0,
        .requiredLimits = null,
        .deviceLostCallback = deviceLostCallback,
        .defaultQueue = .{
            .label = "Default Queue",
            .nextInChain = null,
        },
    };
    var device_user_data = RequestDeviceCallbackData{};
    // When this function returns we know that our callback function was called (there's no need to await it).
    c.wgpuAdapterRequestDevice(adapter, &device_desc, requestDeviceCallback, &device_user_data);

    if (device_user_data.device == null) {
        return error.RequestWebGPUDevice;
    }

    const device = device_user_data.device.?;

    c.wgpuDeviceSetUncapturedErrorCallback(device, deviceErrorCallback, null);

    return device;
}

fn getBindGroupLayout(device: c.WGPUDevice) c.WGPUBindGroupLayout {
    // Consts buffer.
    const binding_0 = c.WGPUBindGroupLayoutEntry{
        .binding = 0,
        .visibility = c.WGPUShaderStage_Compute,
        .buffer = .{
            .type = c.WGPUBufferBindingType_Uniform,
        },
    };

    // Input buffer.
    const binding_1 = c.WGPUBindGroupLayoutEntry{
        .binding = 1,
        .visibility = c.WGPUShaderStage_Compute,
        .buffer = .{
            .type = c.WGPUBufferBindingType_ReadOnlyStorage,
        },
    };

    // Output buffer.
    const binding_2 = c.WGPUBindGroupLayoutEntry{
        .binding = 2,
        .visibility = c.WGPUShaderStage_Compute,
        .buffer = .{
            .type = c.WGPUBufferBindingType_Storage,
        },
    };

    const bindings = [_]c.WGPUBindGroupLayoutEntry{ binding_0, binding_1, binding_2 };

    const layout_desc = c.WGPUBindGroupLayoutDescriptor{
        .label = "Bind Group Layout",
        .nextInChain = null,
        .entryCount = bindings.len,
        .entries = &bindings,
    };

    return c.wgpuDeviceCreateBindGroupLayout(device, &layout_desc);
}

fn getComputePipeline(device: c.WGPUDevice, bind_group_layout: c.WGPUBindGroupLayout) struct { layout: c.WGPUPipelineLayout, pipeline: c.WGPUComputePipeline } {
    const size = 64; // TODO: Get form the outside.
    const n = std.fmt.comptimePrint("{}", .{size});

    const shader_code = @embedFile("./shader.wgsl");

    // Make buffer twice as big to ensure that it's big enough.
    var buffer = [_]u8{0} ** (shader_code.len * 2);
    _ = std.mem.replace(u8, shader_code, "{n}", n, &buffer);

    // Ensure that buffer which hold the code as a string is null-terminated (has a 0 as the last byte).
    const code = buffer[0 .. buffer.len - 1 :0];

    // Load compute shader.
    const shader_code_desc = c.WGPUShaderModuleWGSLDescriptor{
        .chain = .{ .next = null, .sType = c.WGPUSType_ShaderModuleWGSLDescriptor },
        .code = code,
    };

    const shader_module_desc = c.WGPUShaderModuleDescriptor{
        .label = "Shader Module",
        .nextInChain = @ptrCast(&shader_code_desc.chain),
    };

    const compute_shader_module = c.wgpuDeviceCreateShaderModule(device, &shader_module_desc);

    // Create compute pipeline layout.
    const pipeline_layout_desc = c.WGPUPipelineLayoutDescriptor{
        .label = "Pipeline Layout",
        .bindGroupLayoutCount = 1,
        .bindGroupLayouts = &[_]c.WGPUBindGroupLayout{bind_group_layout},
    };
    const layout = c.wgpuDeviceCreatePipelineLayout(device, &pipeline_layout_desc);

    // Create compute pipeline.
    const compute_pipeline_desc = c.WGPUComputePipelineDescriptor{
        .label = "Compute Pipeline",
        .compute = .{
            .constantCount = 0,
            .constants = null,
            .entryPoint = "main",
            .module = compute_shader_module,
        },
        .layout = layout,
    };
    const pipeline = c.wgpuDeviceCreateComputePipeline(device, &compute_pipeline_desc);

    return .{
        .layout = layout,
        .pipeline = pipeline,
    };
}

fn getBuffers(device: c.WGPUDevice, size: u32) struct { consts: c.WGPUBuffer, input: c.WGPUBuffer, output: c.WGPUBuffer, map: c.WGPUBuffer } {
    // Consts buffer.
    const consts_buffer_desc = c.WGPUBufferDescriptor{
        .label = "Consts Buffer",
        .nextInChain = null,
        .mappedAtCreation = 0,
        .size = size,
        .usage = c.WGPUBufferUsage_Uniform | c.WGPUBufferUsage_CopyDst,
    };
    const consts = c.wgpuDeviceCreateBuffer(device, &consts_buffer_desc);

    // Input buffer.
    const input_buffer_desc = c.WGPUBufferDescriptor{
        .label = "Input Buffer",
        .nextInChain = null,
        .mappedAtCreation = 0,
        .size = size,
        .usage = c.WGPUBufferUsage_Storage | c.WGPUBufferUsage_CopyDst,
    };
    const input = c.wgpuDeviceCreateBuffer(device, &input_buffer_desc);

    // Output buffer.
    const output_buffer_desc = c.WGPUBufferDescriptor{
        .label = "Output Buffer",
        .nextInChain = null,
        .mappedAtCreation = 0,
        .size = size,
        .usage = c.WGPUBufferUsage_Storage | c.WGPUBufferUsage_CopySrc,
    };
    const output = c.wgpuDeviceCreateBuffer(device, &output_buffer_desc);

    // Map buffer.
    // This is the intermediate buffer to which we copy the output. This buffer
    // can then be used for reading into the CPU's memory.
    const map_buffer_desc = c.WGPUBufferDescriptor{
        .label = "Map Buffer",
        .nextInChain = null,
        .mappedAtCreation = 0,
        .size = size,
        .usage = c.WGPUBufferUsage_CopyDst | c.WGPUBufferUsage_MapRead,
    };
    const map = c.wgpuDeviceCreateBuffer(device, &map_buffer_desc);

    return .{
        .consts = consts,
        .input = input,
        .output = output,
        .map = map,
    };
}

fn getBindGroup(device: c.WGPUDevice, bind_group_layout: c.WGPUBindGroupLayout, buffer_size: u32, consts_buffer: c.WGPUBuffer, input_buffer: c.WGPUBuffer, output_buffer: c.WGPUBuffer) c.WGPUBindGroup {
    // Consts buffer.
    const entry_0 = c.WGPUBindGroupEntry{
        .binding = 0,
        .buffer = consts_buffer,
        .offset = 0,
        .size = buffer_size,
    };

    // Input buffer.
    const entry_1 = c.WGPUBindGroupEntry{
        .binding = 1,
        .buffer = input_buffer,
        .offset = 0,
        .size = buffer_size,
    };

    // Output buffer.
    const entry_2 = c.WGPUBindGroupEntry{
        .binding = 2,
        .buffer = output_buffer,
        .offset = 0,
        .size = buffer_size,
    };

    const entries = [_]c.WGPUBindGroupEntry{ entry_0, entry_1, entry_2 };

    const bind_group_desc = c.WGPUBindGroupDescriptor{
        .label = "Bind Group",
        .layout = bind_group_layout,
        .entryCount = entries.len,
        .entries = &entries,
    };

    return c.wgpuDeviceCreateBindGroup(device, &bind_group_desc);
}

// Adapter.
const RequestAdapterCallbackData = struct { adapter: ?c.WGPUAdapter = null };

fn requestAdapterCallback(status: c.WGPURequestAdapterStatus, adapter: c.WGPUAdapter, message: [*c]const u8, user_data: ?*anyopaque) callconv(.C) void {
    _ = message;
    var data: *RequestAdapterCallbackData = @ptrCast(@alignCast(user_data));

    if (status == c.WGPURequestAdapterStatus_Success) {
        data.adapter = adapter;
    }
}

// Device.
const RequestDeviceCallbackData = struct { device: ?c.WGPUDevice = null };

fn deviceLostCallback(reason: c.WGPUDeviceLostReason, message: [*c]const u8, user_data: ?*anyopaque) callconv(.C) void {
    _ = user_data;
    const msg = std.mem.span(message);
    std.log.err("Device error: reason {d}\n", .{reason});

    if (msg.len != 0) {
        std.log.err("\tmessage: {s}\n", .{msg});
    }
}

fn deviceErrorCallback(error_type: c.WGPUErrorType, message: [*c]const u8, user_data: ?*anyopaque) callconv(.C) void {
    _ = user_data;
    const msg = std.mem.span(message);
    std.log.err("Device error: type {d}\n", .{error_type});

    if (msg.len != 0) {
        std.log.err("\tmessage: {s}\n", .{msg});
    }
}

fn requestDeviceCallback(status: c.WGPURequestDeviceStatus, device: c.WGPUDevice, message: [*c]const u8, user_data: ?*anyopaque) callconv(.C) void {
    _ = message;
    var data: *RequestDeviceCallbackData = @ptrCast(@alignCast(user_data));

    if (status == c.WGPURequestDeviceStatus_Success) {
        data.device = device;
    }
}

// Buffer.
const BufferMapCallbackData = struct { buffer: c.WGPUBuffer, buffer_size: u32, output: ?f32 = null };

fn bufferMapCallback(status: c.WGPUBufferMapAsyncStatus, user_data: ?*anyopaque) callconv(.C) void {
    const data: *BufferMapCallbackData = @ptrCast(@alignCast(user_data));

    if (status == c.WGPUBufferMapAsyncStatus_Success) {
        const mapped_range_ptr = c.wgpuBufferGetConstMappedRange(data.buffer, 0, data.buffer_size);
        defer c.wgpuBufferUnmap(data.buffer);

        const output: [*]const f32 = @ptrCast(@alignCast(mapped_range_ptr));

        std.debug.print("{s}", .{"["});

        for (0..(data.buffer_size / @sizeOf(f32))) |i| {
            std.debug.print("{d}", .{output[i]});
            if (i < (data.buffer_size / @sizeOf(f32)) - 1) {
                std.debug.print("{s}", .{", "});
            }
        }

        std.debug.print("{s}", .{"]"});

        data.output = output[0];
    }
}

// test "Application - e2e" {
//     const buffer_size = 64 * @sizeOf(f32);

//     const app = try Application(buffer_size).init();
//     defer app.deinit();

//     app.run();

//     // std.debug.print("{}\n", .{app});
// }
