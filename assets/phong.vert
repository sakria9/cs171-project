#version 430 core

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 uv;

out vec3 world_position_frag;
out vec3 normal_frag;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main() {
  world_position_frag = vec3(model * vec4(position, 1));
  normal_frag = mat3(transpose(inverse(model))) * normal;
  gl_Position = projection * view * vec4(world_position_frag, 1.0);
}
