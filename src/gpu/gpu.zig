const std = @import("std");
const testing = std.testing;
const c = @cImport({
    @cInclude("webgpu/wgpu.h");
});

pub fn GPU(comptime n: i32, comptime t: type) type {
    return struct {
        const Self = @This();

        device: c.WGPUDevice = null,
        queue: c.WGPUQueue = null,
        bind_group: c.WGPUBindGroup = null,
        bind_group_layout: c.WGPUBindGroupLayout = null,
        code: [*:0]const u8 = "",
        map_buffer: c.WGPUBuffer = null,
        twiddles_buffer: c.WGPUBuffer = null,
        input_a_buffer: c.WGPUBuffer = null,
        input_b_buffer: c.WGPUBuffer = null,
        output_buffer: c.WGPUBuffer = null,
        buffer_size: u32 = n * @sizeOf(t),

        pub fn init(shader_code: [*:0]const u8) Self {
            var self = Self{};

            self.setDevice();
            self.setQueue();
            self.setBindGroupLayout();
            self.setBuffers();
            self.setBindGroup();

            self.code = shader_code;

            return self;
        }

        pub fn deinit(self: Self) void {
            c.wgpuBindGroupRelease(self.bind_group);
            c.wgpuBufferDestroy(self.map_buffer);
            c.wgpuBufferRelease(self.map_buffer);
            c.wgpuBufferDestroy(self.twiddles_buffer);
            c.wgpuBufferRelease(self.twiddles_buffer);
            c.wgpuBufferDestroy(self.input_a_buffer);
            c.wgpuBufferRelease(self.input_a_buffer);
            c.wgpuBufferDestroy(self.input_b_buffer);
            c.wgpuBufferRelease(self.input_b_buffer);
            c.wgpuBufferDestroy(self.output_buffer);
            c.wgpuBufferRelease(self.output_buffer);
            c.wgpuBindGroupLayoutRelease(self.bind_group_layout);
            c.wgpuQueueRelease(self.queue);
            c.wgpuDeviceRelease(self.device);
        }

        // TODO: Change signature so result and callback are at the front.
        pub fn run(self: Self, func: [*:0]const u8, twiddles: []const t, input_a: []const t, input_b: []const t, result: []t, callback: Callback) !void {
            if (!self.isInitialized()) return error.NotInitialized;

            // Length of inputs and twiddles must equal n.
            if (twiddles.len != input_a.len or input_a.len != input_b.len or input_b.len != n) {
                return error.InvalidLength;
            }

            // Load compute shader.
            const shader_code_desc = c.WGPUShaderModuleWGSLDescriptor{
                .chain = .{
                    .next = null,
                    .sType = c.WGPUSType_ShaderModuleWGSLDescriptor,
                },
                .code = self.code,
            };

            const shader_module_desc = c.WGPUShaderModuleDescriptor{
                .label = "Shader Module",
                .nextInChain = @ptrCast(&shader_code_desc.chain),
            };

            const compute_shader_module = c.wgpuDeviceCreateShaderModule(self.device, &shader_module_desc);

            // Create compute pipeline layout.
            const pipeline_layout_desc = c.WGPUPipelineLayoutDescriptor{
                .label = "Pipeline Layout",
                .bindGroupLayoutCount = 1,
                .bindGroupLayouts = &[_]c.WGPUBindGroupLayout{self.bind_group_layout},
            };
            const layout = c.wgpuDeviceCreatePipelineLayout(self.device, &pipeline_layout_desc);
            defer c.wgpuPipelineLayoutRelease(layout);

            // Create compute pipeline.
            const compute_pipeline_desc = c.WGPUComputePipelineDescriptor{
                .label = "Compute Pipeline",
                .nextInChain = null,
                .compute = .{
                    .constantCount = 0,
                    .constants = null,
                    .entryPoint = func,
                    .module = compute_shader_module,
                },
                .layout = layout,
            };
            const pipeline = c.wgpuDeviceCreateComputePipeline(self.device, &compute_pipeline_desc);
            defer c.wgpuComputePipelineRelease(pipeline);

            // Write into buffers.
            c.wgpuQueueWriteBuffer(self.queue, self.twiddles_buffer, 0, @ptrCast(twiddles), self.buffer_size);
            c.wgpuQueueWriteBuffer(self.queue, self.input_a_buffer, 0, @ptrCast(input_a), self.buffer_size);
            c.wgpuQueueWriteBuffer(self.queue, self.input_b_buffer, 0, @ptrCast(input_b), self.buffer_size);

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
            c.wgpuComputePassEncoderSetPipeline(compute_pass, pipeline);
            c.wgpuComputePassEncoderSetBindGroup(compute_pass, 0, self.bind_group, 0, null);

            // TODO: Update.
            const invocation_count = self.buffer_size / @sizeOf(t);
            const workgroup_size = 32;
            // ceil(invocation_count / workgroup_size)
            const workgroup_count = ((invocation_count + workgroup_size) - 1) / workgroup_size;
            c.wgpuComputePassEncoderDispatchWorkgroups(compute_pass, workgroup_count, 1, 2);

            // Finalize compute pass.
            c.wgpuComputePassEncoderEnd(compute_pass);

            // Before encoder finished.
            c.wgpuCommandEncoderCopyBufferToBuffer(command_encoder, self.output_buffer, 0, self.map_buffer, 0, self.buffer_size);

            // Encode and submit commands to GPU.
            const command_buffer_desc = c.WGPUCommandBufferDescriptor{
                .label = "Command Buffer",
                .nextInChain = null,
            };
            const commands = c.wgpuCommandEncoderFinish(command_encoder, &command_buffer_desc);
            defer c.wgpuCommandBufferRelease(commands);

            c.wgpuQueueSubmit(self.queue, 1, &commands);

            // Read back result.
            var buffer_user_data = BufferMapCallbackData{
                .buffer = self.map_buffer,
                .buffer_size = self.buffer_size,
                .result = @ptrCast(result),
                .callback = callback,
            };

            c.wgpuBufferMapAsync(self.map_buffer, c.WGPUMapMode_Read, 0, self.buffer_size, bufferMapCallback, &buffer_user_data);

            while (!buffer_user_data.is_done) {
                // Mimmicking the `processEvent` function as wgpu-native doesn't come with it yet.
                c.wgpuQueueSubmit(self.queue, 0, null);
            }
        }

        fn setDevice(self: *Self) void {
            // Create WebGPU instance.
            const instance_desc = c.WGPUInstanceDescriptor{
                .nextInChain = null,
            };
            const instance = c.wgpuCreateInstance(&instance_desc) orelse unreachable;
            // The instance isn't needed once the adapter is selected.
            // See: https://eliemichel.github.io/LearnWebGPU/getting-started/adapter-and-device/the-adapter.html#destruction
            defer c.wgpuInstanceRelease(instance);

            // Create WebGPU adapter.
            const adapter_opts = c.WGPURequestAdapterOptions{
                .nextInChain = null,
            };
            var adapter_user_data = RequestAdapterCallbackData{};
            // When this function returns we know that our callback function was called (there's no need to await it).
            c.wgpuInstanceRequestAdapter(instance, &adapter_opts, requestAdapterCallback, &adapter_user_data);

            // The adapter isn't needed once the device is selected.
            // See: https://eliemichel.github.io/LearnWebGPU/getting-started/adapter-and-device/the-device.html#device-request
            const adapter = adapter_user_data.adapter orelse unreachable;
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

            const device = device_user_data.device orelse unreachable;

            c.wgpuDeviceSetUncapturedErrorCallback(device, deviceErrorCallback, null);

            self.device = device;
        }

        fn setQueue(self: *Self) void {
            self.queue = c.wgpuDeviceGetQueue(self.device);
        }

        fn setBindGroupLayout(self: *Self) void {
            // Twiddles buffer.
            const entry_0 = c.WGPUBindGroupLayoutEntry{
                .binding = 0,
                .visibility = c.WGPUShaderStage_Compute,
                .buffer = .{
                    .type = c.WGPUBufferBindingType_Uniform,
                },
            };

            // Input A Buffer.
            const entry_1 = c.WGPUBindGroupLayoutEntry{
                .binding = 1,
                .visibility = c.WGPUShaderStage_Compute,
                .buffer = .{
                    .type = c.WGPUBufferBindingType_ReadOnlyStorage,
                },
            };

            // Input B Buffer.
            const entry_2 = c.WGPUBindGroupLayoutEntry{
                .binding = 2,
                .visibility = c.WGPUShaderStage_Compute,
                .buffer = .{
                    .type = c.WGPUBufferBindingType_ReadOnlyStorage,
                },
            };

            // Output buffer.
            const entry_3 = c.WGPUBindGroupLayoutEntry{
                .binding = 3,
                .visibility = c.WGPUShaderStage_Compute,
                .buffer = .{
                    .type = c.WGPUBufferBindingType_Storage,
                },
            };

            const entries = [_]c.WGPUBindGroupLayoutEntry{ entry_0, entry_1, entry_2, entry_3 };

            const layout_desc = c.WGPUBindGroupLayoutDescriptor{
                .label = "Bind Group Layout",
                .nextInChain = null,
                .entryCount = entries.len,
                .entries = &entries,
            };

            self.bind_group_layout = c.wgpuDeviceCreateBindGroupLayout(self.device, &layout_desc);
        }

        fn setBuffers(self: *Self) void {
            // Twiddles buffer.
            const twiddles_buffer_desc = c.WGPUBufferDescriptor{
                .label = "Twiddles Buffer",
                .nextInChain = null,
                .mappedAtCreation = 0,
                .size = self.buffer_size,
                .usage = c.WGPUBufferUsage_Uniform | c.WGPUBufferUsage_CopyDst,
            };
            self.twiddles_buffer = c.wgpuDeviceCreateBuffer(self.device, &twiddles_buffer_desc);

            // Input A buffer.
            const input_a_buffer_desc = c.WGPUBufferDescriptor{
                .label = "Input A Buffer",
                .nextInChain = null,
                .mappedAtCreation = 0,
                .size = self.buffer_size,
                .usage = c.WGPUBufferUsage_Storage | c.WGPUBufferUsage_CopyDst,
            };
            self.input_a_buffer = c.wgpuDeviceCreateBuffer(self.device, &input_a_buffer_desc);

            // Input B buffer.
            const input_b_buffer_desc = c.WGPUBufferDescriptor{
                .label = "Input B Buffer",
                .nextInChain = null,
                .mappedAtCreation = 0,
                .size = self.buffer_size,
                .usage = c.WGPUBufferUsage_Storage | c.WGPUBufferUsage_CopyDst,
            };
            self.input_b_buffer = c.wgpuDeviceCreateBuffer(self.device, &input_b_buffer_desc);

            // Output buffer.
            const output_buffer_desc = c.WGPUBufferDescriptor{
                .label = "Output Buffer",
                .nextInChain = null,
                .mappedAtCreation = 0,
                .size = self.buffer_size,
                .usage = c.WGPUBufferUsage_Storage | c.WGPUBufferUsage_CopySrc,
            };
            self.output_buffer = c.wgpuDeviceCreateBuffer(self.device, &output_buffer_desc);

            // Map buffer.
            const map_buffer_desc = c.WGPUBufferDescriptor{
                .label = "Map Buffer",
                .nextInChain = null,
                .mappedAtCreation = 0,
                .size = self.buffer_size,
                .usage = c.WGPUBufferUsage_CopyDst | c.WGPUBufferUsage_MapRead,
            };
            self.map_buffer = c.wgpuDeviceCreateBuffer(self.device, &map_buffer_desc);
        }

        fn setBindGroup(self: *Self) void {
            // Twiddles buffer.
            const entry_0 = c.WGPUBindGroupEntry{
                .binding = 0,
                .buffer = self.twiddles_buffer,
                .offset = 0,
                .size = self.buffer_size,
            };

            // Input A buffer.
            const entry_1 = c.WGPUBindGroupEntry{
                .binding = 1,
                .buffer = self.input_a_buffer,
                .offset = 0,
                .size = self.buffer_size,
            };

            // Input B buffer.
            const entry_2 = c.WGPUBindGroupEntry{
                .binding = 2,
                .buffer = self.input_b_buffer,
                .offset = 0,
                .size = self.buffer_size,
            };

            // Output buffer.
            const entry_3 = c.WGPUBindGroupEntry{
                .binding = 3,
                .buffer = self.output_buffer,
                .offset = 0,
                .size = self.buffer_size,
            };

            const entries = [_]c.WGPUBindGroupEntry{ entry_0, entry_1, entry_2, entry_3 };

            const group_desc = c.WGPUBindGroupDescriptor{
                .nextInChain = null,
                .label = "Bind Group",
                .layout = self.bind_group_layout,
                .entryCount = entries.len,
                .entries = &entries,
            };

            self.bind_group = c.wgpuDeviceCreateBindGroup(self.device, &group_desc);
        }

        fn isInitialized(self: Self) bool {
            return self.device != null and
                self.queue != null and
                self.bind_group != null and
                self.bind_group_layout != null and
                std.mem.len(self.code) != 0 and
                self.map_buffer != null and
                self.twiddles_buffer != null and
                self.input_a_buffer != null and
                self.input_b_buffer != null and
                self.output_buffer != null and
                self.buffer_size != 0;
        }
    };
}

const RequestAdapterCallbackData = struct {
    adapter: ?c.WGPUAdapter = null,
};

fn requestAdapterCallback(status: c.WGPURequestAdapterStatus, adapter: c.WGPUAdapter, message: [*c]const u8, user_data: ?*anyopaque) callconv(.C) void {
    _ = message;
    var data: *RequestAdapterCallbackData = @ptrCast(@alignCast(user_data));

    if (status == c.WGPURequestAdapterStatus_Success) {
        data.adapter = adapter;
    }
}

const RequestDeviceCallbackData = struct {
    device: ?c.WGPUDevice = null,
};

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

const BufferMapCallbackData = struct {
    buffer: c.WGPUBuffer,
    buffer_size: u32,
    callback: Callback,
    result: ResultSlice,
    is_done: bool = false,
};

fn bufferMapCallback(status: c.WGPUBufferMapAsyncStatus, user_data: ?*anyopaque) callconv(.C) void {
    const data: *BufferMapCallbackData = @ptrCast(@alignCast(user_data));

    if (status == c.WGPUBufferMapAsyncStatus_Success) {
        const mapped_range_ptr = c.wgpuBufferGetConstMappedRange(data.buffer, 0, data.buffer_size);
        defer c.wgpuBufferUnmap(data.buffer);

        const ctx = Context{
            .buffer_size = data.buffer_size,
            .mapped_range_ptr = mapped_range_ptr,
        };

        cbFunc(ctx, data.result);
    }

    data.is_done = true;
}

pub const Context = struct {
    buffer_size: u32,
    mapped_range_ptr: MappedRangePtr,
};

pub const ResultSlice = *anyopaque;

const MappedRangePtr = ?*const anyopaque;

const Callback = *const fn (Context, ResultSlice) void;

// This function needs to be implemented for NTT and FFT.
fn cbFunc(ctx: Context, result_slice: ResultSlice) void {
    const t = i32;

    const result: [*]t = @ptrCast(@alignCast(result_slice));
    const output: [*]const t = @ptrCast(@alignCast(ctx.mapped_range_ptr));

    for (0..(ctx.buffer_size / @sizeOf(t))) |i| {
        result[i] = output[i];
    }
}

pub fn populateShaderCode(
    comptime shader_code: [*:0]const u8,
    comptime n: i32,
    comptime n_inverse: i32,
    output_buffer: []u8,
) [*:0]const u8 {
    const needle_n = "{n}";
    const needle_n_inverse = "{n_inverse}";
    const replacement_n = std.fmt.comptimePrint("{}", .{n});
    const replacement_n_inverse = std.fmt.comptimePrint("{}", .{n_inverse});

    _ = std.mem.replace(u8, std.mem.span(shader_code), needle_n, replacement_n, output_buffer);
    _ = std.mem.replace(u8, output_buffer, needle_n_inverse, replacement_n_inverse, output_buffer);

    // TODO: Inject a 0 at the end of the output_buffer to ensure that next checks pass.

    const end = std.mem.indexOf(u8, output_buffer, &[_]u8{0}) orelse output_buffer.len;
    const result: [*:0]const u8 = output_buffer[0..end :0];

    return result;
}

test "gpu" {
    const n = 4;
    const n_inverse = 12;
    const t = i32;

    const embedded_shader_code = @embedFile("../gpu.wgsl");
    // const shader_code: [*:0]const u8 = @embedFile("../ntt.wgsl");
    const shader_code: [*:0]const u8 = "const foo : i32 = 42; \n\n" ++ embedded_shader_code;

    // TODO: This is rather wasteful and can cause stack overflows.
    // Use `std.mem.replacementSize` instead to determine the correct size of the buffer.
    var code = [_]u8{0} ** (std.mem.len(shader_code) + 500);
    const populated = populateShaderCode(shader_code, n, n_inverse, &code);

    const gpu = GPU(n, t).init(populated);
    defer gpu.deinit();

    // ---

    const twiddles = [_]t{ 0, 0, 0, 0 };
    const input_a = [_]t{ 1, 2, 3, 4 };
    const input_b = [_]t{ 0, 0, 0, 0 };

    const expected = [_]t{ 1, 2, 3, 4 };
    var result = [_]t{0} ** n;
    try gpu.run("add_twiddles", &twiddles, &input_a, &input_b, &result, cbFunc);

    std.debug.print("{any}", .{result});

    try testing.expectEqualSlices(t, &expected, &result);
}
