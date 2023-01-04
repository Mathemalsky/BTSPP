#include "draw/gui.hpp"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

// imgui library
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "draw/definitions.hpp"
#include "draw/variables.hpp"

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

void initImGuiWindows() {
  // settings window
  imguiwindow::SHOW_SETTINGS_WINDOW = imguiwindow::INITIAL_SHOW_SETTINGS_WINDOW;
  for (const ProblemType& type : problemType::PROBLEM_TYPES) {
    imguiwindow::ACTIVE[type]    = imguiwindow::INITIAL_ACTIVENESS.at(type);
    imguiwindow::COLOR[type]     = imguiwindow::INITIAL_COLOR.at(type);
    imguiwindow::THICKNESS[type] = imguiwindow::INITIAL_THICKNESS.at(type);
  }
  imguiwindow::VERTEX_COLOR = imguiwindow::INITIAL_VERTEX_COLOR;
}

void drawImgui() {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  if (imguiwindow::SHOW_SETTINGS_WINDOW) {
    ImGui::Begin("Settings", &imguiwindow::SHOW_SETTINGS_WINDOW);
    ImGui::Text("mouse x = %f", input::mouse::x);  // DEBUG
    ImGui::Text("mouse y = %f", input::mouse::y);  // DEBUG
    ImGui::Text(
      "Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

    ImGui::Spacing();
    ImGui::Checkbox("BTSP exact", &imguiwindow::ACTIVE[ProblemType::BTSP_exact]);
    ImGui::SliderFloat("thickness##BTSP exact", &imguiwindow::THICKNESS[ProblemType::BTSP_exact], 0.0f, 20.0f, "%.1f");
    ImGui::ColorEdit3("##BTSP exact", (float*) &imguiwindow::COLOR[ProblemType::BTSP_exact]);

    ImGui::Checkbox("TSP  exact", &imguiwindow::ACTIVE[ProblemType::TSP_exact]);
    ImGui::ColorEdit3("##TSP exact", (float*) &imguiwindow::COLOR[ProblemType::TSP_exact]);
    ImGui::SliderFloat("thickness##TSP exact", &imguiwindow::THICKNESS[ProblemType::TSP_exact], 0.0f, 30.0f);

    ImGui::ColorEdit3("vertex color", (float*) &imguiwindow::VERTEX_COLOR);
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
