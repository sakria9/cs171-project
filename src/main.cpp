#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "cloth.h"
#include "input.h"
#include "particle.h"
#include "scene.h"
#include <cstdint>
#include <memory>
#include <stb_image_write.h>

auto drop_center_2d() {
  const Float L = .5f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(0, L, L), 2);
  return particle_system;
}

auto drop_left_2d() {
  const Float L = 1.0f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(L, 2 * L, L), 10);
  return particle_system;
}

auto drop_center_3d() {
  const Float L = .3f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(0, L, L), 2);
  return particle_system;
}

auto drop_center_3d_large() {
  const Float L = 1.0f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(0, L, L), 3);
  return particle_system;
}

auto drop_left_3d() {
  const Float L = 0.3f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(L, 2 * L, L), 3);
  return particle_system;
}

auto drop_left_3d_large() {
  // set particle radius to 0.07 may get better result
  const Float L = 1.0f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(L, 2 * L, L), 9);
  return particle_system;
}

auto drop_left_3d_Large() {
  // CPU: 3.5s
  // GPU: 1s
  const Float L = 2.0f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(1.2f * L, 3.0f, L), 20);
  return particle_system;
}

auto drop_3d_Large() {
  const Float L = 2.5f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(L, 1.0f, 2 * L), 15);
  particle_system->generateParticles(Vec3(-1.5f * L, 4.0f, L), 30);
  return particle_system;
}

int main(int argc, char *argv[]) {
  auto particle_system = drop_left_3d();
#ifndef USE_CPU
  particle_system->external_pcisph_init(); // comment this line to use CPU
  std::cerr << "Using GPU" << std::endl;
#else
  std::cerr << "Using CPU" << std::endl;
#endif
  if (argc == 2) {
    // read from file
    std::string filepath = argv[1];
    std::cerr << "Read from " << filepath << std::endl;
    std::fstream file;
    file.open(filepath, std::ios::in | std::ios::binary);
    size_t n;
    file.read((char *)&n, sizeof(int));
    size_t rounds;
    file.read((char *)&rounds, sizeof(int));
    std::cerr << "Read " << n << " particles, " << rounds << " frames"
              << std::endl;
    std::vector<Float> data;
    data.resize(3 * n * rounds);
    file.read((char *)data.data(), sizeof(Float) * 3 * n * rounds);
    file.close();
    std::cerr << "Read data " << data.size() << std::endl;
    particle_system->use_data_init(n, std::move(data));
  }
  if (argc == 3) {
    // simulate and output to file
    std::string filepath = argv[1];
    size_t rounds = std::stoull(argv[2]);

    std::vector<Float> data;
    size_t n = particle_system->particles.size();
    data.resize(3 * n);

    std::cerr << "Output to " << filepath << std::endl;
    std::fstream file;
    file.open(filepath, std::ios::out | std::ios::binary);
    file.write((char *)&n, sizeof(int));
    file.write((char *)&rounds, sizeof(int));

    for (int cnt = 0; cnt < rounds; cnt++) {
      auto st = std::chrono::high_resolution_clock::now();
      std::cerr << "rendering frame " << cnt << std::endl;
      particle_system->fixedUpdate();
      for (int i = 0; i < n; i++) {
        data[3 * i + 0] = particle_system->particles[i].x.x;
        data[3 * i + 1] = particle_system->particles[i].x.y;
        data[3 * i + 2] = particle_system->particles[i].x.z;
      }
      file.write((char *)data.data(), sizeof(Float) * 3 * n);
      auto ed = std::chrono::high_resolution_clock::now();
      auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(ed - st);
      std::cerr << "frame " << cnt << " took " << dur.count() << "ms"
                << std::endl;
    }
    file.close();
    exit(0);
  }

  // window
  constexpr int window_width = 1920;
  constexpr int window_height = 1080;

  // cloth
  constexpr Float cloth_weight = Float(2);
  constexpr UVec2 mass_dim = {40, 30};
  constexpr Float dx_local = Float(0.1);
  constexpr Float stiffness = Float(15);
  constexpr Float damping_ratio = Float(0.0015);
  std::vector<IVec2> fixed_masses{{0, -1}, {-1, -1}};

  /// setup window
  GLFWwindow *window;
  {
    if (!glfwInit()) // initialize glfw library
      return -1;

    // setting glfw window hints and global configurations
    {
      glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
      glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
      glfwWindowHint(GLFW_OPENGL_PROFILE,
                     GLFW_OPENGL_CORE_PROFILE); // use core mode
      // glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE); // use debug
      // context
#ifdef __APPLE__
      glfwWindowHint(
          GLFW_OPENGL_FORWARD_COMPAT,
          GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif
    }

    // create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(window_width, window_height,
                              "CS171 Final Project", NULL, NULL);
    if (!window) {
      glfwTerminate();
      return -1;
    }

    // make the window's context current
    glfwMakeContextCurrent(window);

    // load Opengl
    if (!gladLoadGL()) {
      glfwTerminate();
      return -1;
    }

    // setup call back functions
    glfwSetFramebufferSizeCallback(window, Input::CallBackResizeFlareBuffer);
  }

  /// main Loop
  {
    // shader
    Shader::Initialize();

    // scene
    Scene scene(45);
    scene.camera.transform.position = {0, 0.8, -6};
    scene.camera.transform.rotation = {0, 0, 1, 0};
    scene.light_position = {0, 3, -10};
    scene.light_color = Vec3(1, 1, 1) * Float(1.125);

    particle_system->mesh_sphere =
        std::make_shared<Mesh>(MeshPrimitiveType::sphere);
    {
      auto objs = particle_system->boundryIndicators();
      scene.objects.insert(scene.objects.end(), objs.begin(), objs.end());
    }

    // loop until the user closes the window
    Input::Start(window);
    //    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    uint8_t *pixels = new uint8_t[3 * window_width * window_height];

    std::string prefix = "~/cgout/2d/";
    size_t cnt = 0;
    bool start = false;
    while (!glfwWindowShouldClose(window)) {
      Input::Update();
      Time::Update();
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      static Float last_gen_time = Time::elapsed_time;
      if (Time::elapsed_time - last_gen_time > 0.2f) {
        last_gen_time = Time::elapsed_time;
        // particle_system->generateParticles(0);
      }

      /// terminate
      if (Input::GetKey(KeyCode::Escape))
        glfwSetWindowShouldClose(window, true);

      if (Input::GetKey(KeyCode::Return))
        start = true;

      // /// fixed update
      // for (unsigned i = 0; i < Time::fixed_update_times_this_frame; ++i) {
      if (Input::GetKey(KeyCode::Space) || start) //! only when space is pressed
      {
        scene.FixedUpdate();
        particle_system->fixedUpdate();
      }
      // }

      /// camera update
      { scene.Update(); }

      /// render
      {
        scene.RenderUpdate();
        particle_system->renderParticle(scene);
        // particle_system->renderSurface(scene);
      }

      // swap front and back buffers
      glfwSwapBuffers(window);

      // poll for and process events
      glfwPollEvents();

      // std::cerr << "render time: " << Time::delta_time << std::endl;

      // get window capture
      if (start) {
        // std::cerr << "capture " << cnt << std::endl;
        // glReadPixels(0, 0, window_width, window_height, GL_RGB,
        //              GL_UNSIGNED_BYTE, pixels);
        // std::string filename = prefix + std::to_string(cnt++) + ".png";
        // stbi_flip_vertically_on_write(true);
        // stbi_write_png(filename.c_str(), window_width, window_height, 3,
        // pixels,
        //                window_width * 3);
      }
    }
  }

  glfwTerminate();
  return 0;
}
