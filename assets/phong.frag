#version 430 core

out vec4 frag_color;

in vec3 world_position_frag;
in vec3 normal_frag;

uniform vec3 object_color;
uniform vec3 light_position;
uniform vec3 light_color;
uniform vec3 camera_position;

uniform bool is_bidirectional;

void main() {
  // ambient
  const float AMBIENT_STRENGTH = 0.1;
  vec3 ambient = AMBIENT_STRENGTH * light_color;

  // diffuse
  vec3 normal_frag_normalized = normalize(normal_frag);
  vec3 light_direction = normalize(light_position - world_position_frag);
  vec3 diffuse = is_bidirectional ? abs(dot(normal_frag_normalized, light_direction)) * light_color :
                                    max(dot(normal_frag_normalized, light_direction), 0.0) * light_color;

  // specular
  const float SPECULAR_STRENGTH = 0.5;
  vec3 camera_direction = normalize(camera_position - world_position_frag);
  vec3 reflect_direction = reflect(-light_direction, normal_frag_normalized);
  vec3 specular = is_bidirectional ? SPECULAR_STRENGTH * pow(abs(dot(camera_direction, reflect_direction)), 32) * light_color :
                                     SPECULAR_STRENGTH * pow(max(dot(camera_direction, reflect_direction), 0.0), 32) * light_color;

  vec3 result = (ambient + diffuse + specular) * object_color;
  frag_color = vec4(result, 1.0);
}
