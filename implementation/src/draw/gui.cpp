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
    imguiwindow::ACTIVE[(unsigned int) type]    = imguiwindow::INITIAL_ACTIVENESS.at((unsigned int) type);
    imguiwindow::COLOUR[(unsigned int) type]    = imguiwindow::INITIAL_COLOUR.at((unsigned int) type);
    imguiwindow::THICKNESS[(unsigned int) type] = imguiwindow::INITIAL_THICKNESS.at((unsigned int) type);
  }
  imguiwindow::VERTEX_COLOUR = imguiwindow::INITIAL_VERTEX_COLOUR;
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
    ImGui::Checkbox("BTSP exact", &imguiwindow::ACTIVE[(unsigned int) ProblemType::BTSP_exact]);
    ImGui::SliderFloat(
        "thickness##BTSP exact", &imguiwindow::THICKNESS[(unsigned int) ProblemType::BTSP_exact], 0.0f, 20.0f, "%.1f");
    ImGui::ColorEdit3("##BTSP exact", (float*) &imguiwindow::COLOUR[(unsigned int) ProblemType::BTSP_exact]);

    ImGui::Checkbox("TSP  exact", &imguiwindow::ACTIVE[(unsigned int) ProblemType::TSP_exact]);
    ImGui::ColorEdit3("##TSP exact", (float*) &imguiwindow::COLOUR[(unsigned int) ProblemType::TSP_exact]);
    ImGui::SliderFloat(
        "thickness##TSP exact", &imguiwindow::THICKNESS[(unsigned int) ProblemType::TSP_exact], 0.0f, 30.0f);

    ImGui::ColorEdit3("vertex colour", (float*) &imguiwindow::VERTEX_COLOUR);
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
