#include "mesh.h"

Mesh::Mesh(const MeshPrimitiveType& mesh_primitive) :
    draw_mode(DrawMode::elements),
    buffer_data_usage_vbo(GL_STATIC_DRAW),
    buffer_data_usage_ebo(GL_STATIC_DRAW),
    is_bidirectional(false),
    vao(0),
    vbo(0),
    ebo(0) {
  // setup vertices and indices
  if (mesh_primitive == MeshPrimitiveType::cube) {
    vertices = {
        MeshVertex(Vec3(-0.5, -0.5, -0.5), Vec3( 0.0,  0.0, -1.0), Vec2(0.0, 0.0)), // forward
        MeshVertex(Vec3( 0.5, -0.5, -0.5), Vec3( 0.0,  0.0, -1.0), Vec2(1.0, 0.0)),
        MeshVertex(Vec3( 0.5,  0.5, -0.5), Vec3( 0.0,  0.0, -1.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3( 0.5,  0.5, -0.5), Vec3( 0.0,  0.0, -1.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3(-0.5,  0.5, -0.5), Vec3( 0.0,  0.0, -1.0), Vec2(0.0, 1.0)),
        MeshVertex(Vec3(-0.5, -0.5, -0.5), Vec3( 0.0,  0.0, -1.0), Vec2(0.0, 0.0)),

        MeshVertex(Vec3(-0.5, -0.5,  0.5), Vec3( 0.0,  0.0,  1.0), Vec2(0.0, 0.0)), // back
        MeshVertex(Vec3( 0.5, -0.5,  0.5), Vec3( 0.0,  0.0,  1.0), Vec2(1.0, 0.0)),
        MeshVertex(Vec3( 0.5,  0.5,  0.5), Vec3( 0.0,  0.0,  1.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3( 0.5,  0.5,  0.5), Vec3( 0.0,  0.0,  1.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3(-0.5,  0.5,  0.5), Vec3( 0.0,  0.0,  1.0), Vec2(0.0, 1.0)),
        MeshVertex(Vec3(-0.5, -0.5,  0.5), Vec3( 0.0,  0.0,  1.0), Vec2(0.0, 0.0)),

        MeshVertex(Vec3(-0.5,  0.5,  0.5), Vec3(-1.0,  0.0,  0.0), Vec2(1.0, 0.0)), // left
        MeshVertex(Vec3(-0.5,  0.5, -0.5), Vec3(-1.0,  0.0,  0.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3(-0.5, -0.5, -0.5), Vec3(-1.0,  0.0,  0.0), Vec2(0.0, 1.0)),
        MeshVertex(Vec3(-0.5, -0.5, -0.5), Vec3(-1.0,  0.0,  0.0), Vec2(0.0, 1.0)),
        MeshVertex(Vec3(-0.5, -0.5,  0.5), Vec3(-1.0,  0.0,  0.0), Vec2(0.0, 0.0)),
        MeshVertex(Vec3(-0.5,  0.5,  0.5), Vec3(-1.0,  0.0,  0.0), Vec2(1.0, 0.0)),

        MeshVertex(Vec3( 0.5,  0.5,  0.5), Vec3( 1.0,  0.0,  0.0), Vec2(1.0, 0.0)), // right
        MeshVertex(Vec3( 0.5,  0.5, -0.5), Vec3( 1.0,  0.0,  0.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3( 0.5, -0.5, -0.5), Vec3( 1.0,  0.0,  0.0), Vec2(0.0, 1.0)),
        MeshVertex(Vec3( 0.5, -0.5, -0.5), Vec3( 1.0,  0.0,  0.0), Vec2(0.0, 1.0)),
        MeshVertex(Vec3( 0.5, -0.5,  0.5), Vec3( 1.0,  0.0,  0.0), Vec2(0.0, 0.0)),
        MeshVertex(Vec3( 0.5,  0.5,  0.5), Vec3( 1.0,  0.0,  0.0), Vec2(1.0, 0.0)),

        MeshVertex(Vec3(-0.5, -0.5, -0.5), Vec3( 0.0, -1.0,  0.0), Vec2(0.0, 1.0)), // down
        MeshVertex(Vec3( 0.5, -0.5, -0.5), Vec3( 0.0, -1.0,  0.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3( 0.5, -0.5,  0.5), Vec3( 0.0, -1.0,  0.0), Vec2(1.0, 0.0)),
        MeshVertex(Vec3( 0.5, -0.5,  0.5), Vec3( 0.0, -1.0,  0.0), Vec2(1.0, 0.0)),
        MeshVertex(Vec3(-0.5, -0.5,  0.5), Vec3( 0.0, -1.0,  0.0), Vec2(0.0, 0.0)),
        MeshVertex(Vec3(-0.5, -0.5, -0.5), Vec3( 0.0, -1.0,  0.0), Vec2(0.0, 1.0)),

        MeshVertex(Vec3(-0.5,  0.5, -0.5), Vec3( 0.0,  1.0,  0.0), Vec2(0.0, 1.0)), // up
        MeshVertex(Vec3( 0.5,  0.5, -0.5), Vec3( 0.0,  1.0,  0.0), Vec2(1.0, 1.0)),
        MeshVertex(Vec3( 0.5,  0.5,  0.5), Vec3( 0.0,  1.0,  0.0), Vec2(1.0, 0.0)),
        MeshVertex(Vec3( 0.5,  0.5,  0.5), Vec3( 0.0,  1.0,  0.0), Vec2(1.0, 0.0)),
        MeshVertex(Vec3(-0.5,  0.5,  0.5), Vec3( 0.0,  1.0,  0.0), Vec2(0.0, 0.0)),
        MeshVertex(Vec3(-0.5,  0.5, -0.5), Vec3( 0.0,  1.0,  0.0), Vec2(0.0, 1.0)),
      };
    indices = {
        UVec3( 0,  1,  2), UVec3( 3,  4,  5),
        UVec3( 6,  7,  8), UVec3( 9, 10, 11),
        UVec3(12, 13, 14), UVec3(15, 16, 17),
        UVec3(18, 19, 20), UVec3(21, 22, 23),
        UVec3(24, 25, 26), UVec3(27, 28, 29),
        UVec3(30, 31, 32), UVec3(33, 34, 35),
      };
  }
  else if (mesh_primitive == MeshPrimitiveType::sphere) {
    // compute vertices and indices for +x part
    const Float deg_to_rad = pi / Float(180);
    for (int i = 0; i < primitive_sphere_resolution; ++i) {
      Float phi = deg_to_rad * (Float(45) - Float(90) * Float(i) / (primitive_sphere_resolution - 1)); // [45°, -45°]
      Vec3 d_phi(-Float(std::sin(phi)), Float(std::cos(phi)), 0);
      auto k1 = unsigned(i * primitive_sphere_resolution);
      auto k2 = k1 + unsigned(primitive_sphere_resolution);
      Float t = Float(i) / (primitive_sphere_resolution - 1);
      for (int j = 0; j < primitive_sphere_resolution; ++j, ++k1, ++k2) {
        // vertex
        Float theta = deg_to_rad * (Float(-45) + Float(90) * Float(j) / (primitive_sphere_resolution - 1)); // [-45°, 45°]
        Vec3 d_theta(-Float(std::sin(theta)), 0, -Float(std::cos(theta)));
        Vec3 n = glm::normalize(glm::cross(d_theta, d_phi));
        Float s = (Float)j / (primitive_sphere_resolution - 1);
        vertices.emplace_back(MeshVertex(n * Float(0.5), n, Vec2(s, t)));

        // index
        if (i < (primitive_sphere_resolution - 1) && j < (primitive_sphere_resolution - 1)) {
          indices.emplace_back(UVec3(k1, k2, k1 + 1));
          indices.emplace_back(UVec3(k1 + 1, k2, k2 + 1));
        }
      }
    }
    size_t num_vertices_per_part = vertices.size();
    size_t num_indices_per_part = indices.size();
    // -x
    auto index_start = unsigned(vertices.size());
    for (size_t i = 0; i < num_vertices_per_part; ++i) {
      Vec3 v = vertices[i].position;
      Vec3 n = vertices[i].normal;
      vertices.emplace_back(MeshVertex(Vec3(-v.x, v.y, -v.z), Vec3(-n.x, n.y, -n.z), vertices[i].uv));
    }
    for (size_t i = 0; i < num_indices_per_part; ++i) {
      UVec3 index = indices[i];
      indices.emplace_back(UVec3(index.x + index_start, index.y + index_start, index.z + index_start));
    }
    // +y
    index_start = unsigned(vertices.size());
    for (size_t i = 0; i < num_vertices_per_part; ++i) {
      Vec3 v = vertices[i].position;
      Vec3 n = vertices[i].normal;
      vertices.emplace_back(MeshVertex(Vec3(-v.z, v.x, -v.y), Vec3(-n.z, n.x, -n.y), vertices[i].uv));
    }
    for (size_t i = 0; i < num_indices_per_part; ++i) {
      UVec3 index = indices[i];
      indices.emplace_back(UVec3(index.x + index_start, index.y + index_start, index.z + index_start));
    }
    // -y
    index_start = unsigned(vertices.size());
    for (size_t i = 0; i < num_vertices_per_part; ++i) {
      Vec3 v = vertices[i].position;
      Vec3 n = vertices[i].normal;
      vertices.emplace_back(MeshVertex(Vec3(-v.z, -v.x, v.y), Vec3(-n.z, -n.x, n.y), vertices[i].uv));
    }
    for (size_t i = 0; i < num_indices_per_part; ++i) {
      UVec3 index = indices[i];
      indices.emplace_back(UVec3(index.x + index_start, index.y + index_start, index.z + index_start));
    }
    // +z
    index_start = unsigned(vertices.size());
    for (size_t i = 0; i < num_vertices_per_part; ++i) {
      Vec3 v = vertices[i].position;
      Vec3 n = vertices[i].normal;
      vertices.emplace_back(MeshVertex(Vec3(-v.z, v.y, v.x), Vec3(-n.z, n.y, n.x), vertices[i].uv));
    }
    for (size_t i = 0; i < num_indices_per_part; ++i) {
      UVec3 index = indices[i];
      indices.emplace_back(UVec3(index.x + index_start, index.y + index_start, index.z + index_start));
    }
    // -z
    index_start = unsigned(vertices.size());
    for (size_t i = 0; i < num_vertices_per_part; ++i) {
      Vec3 v = vertices[i].position;
      Vec3 n = vertices[i].normal;
      vertices.emplace_back(MeshVertex(Vec3(v.z, v.y, -v.x), Vec3(n.z, n.y, -n.x), vertices[i].uv));
    }
    for (size_t i = 0; i < num_indices_per_part; ++i) {
      UVec3 index = indices[i];
      indices.emplace_back(UVec3(index.x + index_start, index.y + index_start, index.z + index_start));
    }
  }
  else {
    std::cerr << "Mesh::Mesh(): Invalid mesh primitive" << std::endl;
    abort();
  }

  // construct vao for triangles
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(MeshVertex) * this->vertices.size(), this->vertices.data(), this->buffer_data_usage_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, position));
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, normal));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, uv));
  glEnableVertexAttribArray(2);
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(UVec3) * this->indices.size(), this->indices.data(), this->buffer_data_usage_ebo);
  glBindVertexArray(0);
}

Mesh::Mesh(std::vector<MeshVertex>&& vertices,
           GLenum buffer_data_usage_vbo/* = GL_STATIC_DRAW*/,
           bool is_bidirectional/* = false*/) :
    draw_mode(DrawMode::arrays),
    vertices(std::move(vertices)),
    buffer_data_usage_vbo(buffer_data_usage_vbo),
    buffer_data_usage_ebo(GL_STATIC_DRAW),
    is_bidirectional(is_bidirectional),
    vao(0),
    vbo(0),
    ebo(0) {
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(MeshVertex) * this->vertices.size(), this->vertices.data(), this->buffer_data_usage_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, position));
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, normal));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, uv));
  glEnableVertexAttribArray(2);
  glBindVertexArray(0);
}

Mesh::Mesh(std::vector<MeshVertex>&& vertices, std::vector<UVec3>&& indices,
           GLenum buffer_data_usage_vbo/* = GL_STATIC_DRAW*/, GLenum buffer_data_usage_ebo/* = GL_STATIC_DRAW*/,
           bool is_bidirectional/* = false*/) :
    draw_mode(DrawMode::elements),
    vertices(std::move(vertices)),
    indices(std::move(indices)),
    buffer_data_usage_vbo(buffer_data_usage_vbo),
    buffer_data_usage_ebo(buffer_data_usage_ebo),
    is_bidirectional(is_bidirectional),
    vao(0),
    vbo(0),
    ebo(0) {
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(MeshVertex) * this->vertices.size(), this->vertices.data(), this->buffer_data_usage_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, position));
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, normal));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(MeshVertex), (void*)offsetof(MeshVertex, uv));
  glEnableVertexAttribArray(2);
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(UVec3) * this->indices.size(), this->indices.data(), this->buffer_data_usage_ebo);
  glBindVertexArray(0);
}

void Mesh::SetObject(const std::shared_ptr<Object>& new_object) {
  object = new_object;
}

bool Mesh::IsBidirectional() const {
  return is_bidirectional;
}

void Mesh::DrawTriangles() const {
  glBindVertexArray(vao);
  if (draw_mode == DrawMode::arrays)
    glDrawArrays(GL_TRIANGLES, 0, vertices.size());
  else if (draw_mode == DrawMode::elements)
    glDrawElements(GL_TRIANGLES, 3 * indices.size(), GL_UNSIGNED_INT, 0);
  glBindVertexArray(0);
}

/*virtual*/ void Mesh::FixedUpdate() { }
