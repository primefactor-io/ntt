const std = @import("std");

// Although this function looks imperative, note that its job is to
// declaratively construct a build graph that will be executed by an external
// runner.
pub fn build(b: *std.Build) !void {
    // Standard target options allows the person running `zig build` to choose
    // what target to build for. Here we do not override the defaults, which
    // means any target is allowed, and the default is native. Other options
    // for restricting supported target set are available.
    const target = b.standardTargetOptions(.{});

    // Determine Operating System.
    const os = switch (target.result.os.tag) {
        .linux => "linux",
        .macos => "macos",
        .windows => "windows",
        else => null,
    };

    // Determine CPU Architecture.
    const arch = switch (target.result.cpu.arch) {
        .x86_64 => "x86_64",
        .aarch64 => "aarch64",
        else => null,
    };

    // Exit early if system isn't supported.
    if (os == null or arch == null) {
        return error.UnsupportedSystem;
    }

    // WebGPU paths.
    const wgpuIncludePath = b.fmt("{s}", .{"webgpu/include"});
    const wgpuLibraryPath = b.fmt("webgpu/bin/{s}-{s}", .{ os.?, arch.? });

    // Standard optimization options allow the person running `zig build` to select
    // between Debug, ReleaseSafe, ReleaseFast, and ReleaseSmall. Here we do not
    // set a preferred release mode, allowing the user to decide how to optimize.
    const optimize = b.standardOptimizeOption(.{});

    // This creates a "module", which represents a collection of source files alongside
    // some compilation options, such as optimization mode and linked system libraries.
    // Every executable or library we compile will be based on one or more modules.
    const lib_mod = b.createModule(.{
        // `root_source_file` is the Zig "entry point" of the module. If a module
        // only contains e.g. external object files, you can make this `null`.
        // In this case the main source file is merely a path, however, in more
        // complicated build scripts, this could be a generated file.
        .root_source_file = b.path("src/root.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Now, we will create a static library based on the module we created above.
    // This creates a `std.Build.Step.Compile`, which is the build step responsible
    // for actually invoking the compiler.
    const lib = b.addLibrary(.{
        .linkage = .static,
        .name = "ntt",
        .root_module = lib_mod,
    });

    lib.addIncludePath(b.path(wgpuIncludePath));
    lib.addLibraryPath(b.path(wgpuLibraryPath));
    lib.linkSystemLibrary("wgpu_native");

    // This declares intent for the library to be installed into the standard
    // location when the user invokes the "install" step (the default step when
    // running `zig build`).
    b.installArtifact(lib);

    // Add `check` step for ZLS' Build-On-Save functionality.
    // See: https://zigtools.org/zls/guides/build-on-save/
    const lib_check = b.addLibrary(.{
        .linkage = .static,
        .name = "ntt",
        .root_module = lib_mod,
    });

    lib_check.addIncludePath(b.path(wgpuIncludePath));
    lib_check.addLibraryPath(b.path(wgpuLibraryPath));
    lib_check.linkSystemLibrary("wgpu_native");

    const check = b.step("check", "Check if compilation succeeds");
    check.dependOn(&lib_check.step);

    // Creates a step for unit testing. This only builds the test executable
    // but does not run it.
    const lib_unit_tests = b.addTest(.{
        .root_module = lib_mod,
    });

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);

    // Similar to creating the run step earlier, this exposes a `test` step to
    // the `zig build --help` menu, providing a way for the user to request
    // running the unit tests.
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    // This declares intents for the executables to be installed into the
    // standard location when the user invokes the "install" step (the default
    // step when running `zig build`).
    b.installArtifact(lib_unit_tests);
}
