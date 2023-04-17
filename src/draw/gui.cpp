#include "draw/gui.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// imgui library
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "draw/definitions.hpp"
#include "draw/drawdata.hpp"
#include "draw/variables.hpp"

using namespace drawing;

const char* imguiVersionHints() {
  // GL 4.4 + GLSL 440
  const char* glsl_version = "#version 440";
  glfwWindowHint(GLFW_SAMPLES, 4);  // 4x antialiasing
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 4);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
  return glsl_version;
}

void setUpImgui(GLFWwindow* window, const char* glsl_version) {
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void) io;
  // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls

  // Setup Dear ImGui style
  ImGui::StyleColorsDark();
  // ImGui::StyleColorsClassic();

  // Setup Platform/Renderer backends
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init(glsl_version);
}

void drawImgui(Appearance& appearance) {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  if (drawing::SHOW_DEBUG_WINDOW) {
    ImGui::Begin("Debug", &drawing::SHOW_DEBUG_WINDOW);
    ImGui::Text("Node in motion: %d", input::mouse::NODE_IN_MOTION);
    ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
    ImGui::End();
  }

  if (drawing::SHOW_SETTINGS_WINDOW) {
    ImGui::Begin("Settings", &drawing::SHOW_SETTINGS_WINDOW);

    ImGui::Checkbox("BTSP approx", &drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_approx)]);
    ImGui::ColorEdit3("##BTSP approx", (float*) &appearance.colour[std::to_underlying(ProblemType::BTSP_approx)]);
    ImGui::SliderFloat("thickness##BTSP approx", &appearance.thickness[std::to_underlying(ProblemType::BTSP_approx)], 0.0f, 30.0f, "%.1f");

    ImGui::Checkbox("biconnectd graph##BTSP approx", &drawing::BTSP_DRAW_BICONNECTED_GRAPH);
    ImGui::SameLine();
    ImGui::Checkbox("open ear decomp.", &drawing::BTSP_DRAW_OPEN_EAR_DECOMPOSITION);
    ImGui::Checkbox("hamilton cycle", &drawing::BTSP_DRAW_HAMILTON_CYCLE);

    ImGui::Checkbox("BTSPP approx", &drawing::ACTIVE[std::to_underlying(ProblemType::BTSPP_approx)]);
    ImGui::ColorEdit3("##BTSPP approx", (float*) &appearance.colour[std::to_underlying(ProblemType::BTSPP_approx)]);
    ImGui::SliderFloat("thickness##BTSPP approx", &appearance.thickness[std::to_underlying(ProblemType::BTSPP_approx)], 0.0f, 30.0f,
                       "%.1f");

    ImGui::Checkbox("biconnectd graph##BTSPP approx", &drawing::BTSPP_DRAW_BICONNECTED_GRAPH);
    ImGui::SameLine();
    ImGui::Checkbox("hamilton path", &drawing::BTSPP_DRAW_HAMILTON_PATH);

    ImGui::Checkbox("BTSP exact", &drawing::ACTIVE[std::to_underlying(ProblemType::BTSP_exact)]);
    ImGui::SliderFloat("thickness##BTSP exact", &appearance.thickness[std::to_underlying(ProblemType::BTSP_exact)], 0.0f, 20.0f, "%.1f");
    ImGui::ColorEdit3("##BTSP exact", (float*) &appearance.colour[std::to_underlying(ProblemType::BTSP_exact)]);

    ImGui::Checkbox("BTSPP exact", &drawing::ACTIVE[std::to_underlying(ProblemType::BTSPP_exact)]);
    ImGui::SliderFloat("thickness##BTSPP exact", &appearance.thickness[std::to_underlying(ProblemType::BTSPP_exact)], 0.0f, 20.0f, "%.1f");
    ImGui::ColorEdit3("##BTSPP exact", (float*) &appearance.colour[std::to_underlying(ProblemType::BTSPP_exact)]);

    ImGui::Checkbox("TSP  exact", &drawing::ACTIVE[std::to_underlying(ProblemType::TSP_exact)]);
    ImGui::ColorEdit3("##TSP exact", (float*) &appearance.colour[std::to_underlying(ProblemType::TSP_exact)]);
    ImGui::SliderFloat("thickness##TSP exact", &appearance.thickness[std::to_underlying(ProblemType::TSP_exact)], 0.0f, 30.0f, "%.1f");

    ImGui::ColorEdit4("clear colour", (float*) &appearance.clearColour);
    ImGui::ColorEdit3("vertex colour", (float*) &appearance.vertexColour);
    ImGui::End();
  }

  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void cleanUpImgui() {
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
}
