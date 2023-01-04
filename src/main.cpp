#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "cloth.h"
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
  const Float L = 1.0f;
  auto particle_system = std::make_shared<ParticleSystem>(L);
  particle_system->generateParticles(Vec3(L, 2 * L, L), 9);
  return particle_system;
}

int main() {

  /// settings

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

    auto particle_system = drop_left_3d();
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

      // /// fixed update
      // for (unsigned i = 0; i < Time::fixed_update_times_this_frame; ++i) {
      if (Input::GetKey(KeyCode::Space)) //! only when space is pressed
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
        particle_system->renderSurface(scene);
      }

      // swap front and back buffers
      glfwSwapBuffers(window);

      // poll for and process events
      glfwPollEvents();

      // get window capture
      // std::cerr << "capture " << cnt << std::endl;
      // glReadPixels(0, 0, window_width, window_height, GL_RGB,
      // GL_UNSIGNED_BYTE,
      //              pixels);
      // std::string filename = prefix + std::to_string(cnt++) + ".png";
      // stbi_flip_vertically_on_write(true);
      // stbi_write_png(filename.c_str(), window_width, window_height, 3,
      // pixels,
      //                window_width * 3);
    }
  }

  glfwTerminate();
  return 0;
}
