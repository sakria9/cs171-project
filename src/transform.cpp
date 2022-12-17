#include "transform.h"

Transform::Transform(const Vec3& position, const Quat& rotation, const Vec3& scale) : position(position), rotation(rotation), scale(scale) { }
Transform::Transform() : Transform(Vec3(0, 0, 0), Quat(1, 0, 0, 0), Vec3(1, 1, 1)) { }

Mat3 Transform::RotationMat() const {
  Float
    xx = rotation.x * rotation.x, yy = rotation.y * rotation.y, zz = rotation.z * rotation.z,
    xy = rotation.x * rotation.y, yz = rotation.y * rotation.z, zw = rotation.z * rotation.w, wx = rotation.w * rotation.x,
    xz = rotation.x * rotation.z, yw = rotation.y * rotation.w;
  return {
    1 - 2 * (yy + zz), 2 * (xy + zw), 2 * (xz - yw),
    2 * (xy - zw), 1 - 2 * (xx + zz), 2 * (yz + wx),
    2 * (xz + yw), 2 * (yz - wx), 1 - 2 * (xx + yy) };
}

Mat4 Transform::ModelMat() const {
  Float
    xx = rotation.x * rotation.x, yy = rotation.y * rotation.y, zz = rotation.z * rotation.z,
    xy = rotation.x * rotation.y, yz = rotation.y * rotation.z, zw = rotation.z * rotation.w, wx = rotation.w * rotation.x,
    xz = rotation.x * rotation.z, yw = rotation.y * rotation.w;
  return {
    scale.x * (1 - 2 * (yy + zz)), scale.x * (2 * (xy + zw)), scale.x * (2 * (xz - yw)), 0,
    scale.y * (2 * (xy - zw)), scale.y * (1 - 2 * (xx + zz)), scale.y * (2 * (yz + wx)), 0,
    scale.z * (2 * (xz + yw)), scale.z * (2 * (yz - wx)), scale.z * (1 - 2 * (xx + yy)), 0,
    position.x, position.y, position.z, 1 };
}

Mat4 Transform::ModelMatNotScaled() const {
  Float
    xx = rotation.x * rotation.x, yy = rotation.y * rotation.y, zz = rotation.z * rotation.z,
    xy = rotation.x * rotation.y, yz = rotation.y * rotation.z, zw = rotation.z * rotation.w, wx = rotation.w * rotation.x,
    xz = rotation.x * rotation.z, yw = rotation.y * rotation.w;
  return {
    1 - 2 * (yy + zz), 2 * (xy + zw), 2 * (xz - yw), 0,
    2 * (xy - zw), 1 - 2 * (xx + zz), 2 * (yz + wx), 0,
    2 * (xz + yw), 2 * (yz - wx), 1 - 2 * (xx + yy), 0,
    position.x, position.y, position.z, 1 };
}

Vec3 Transform::Right()   const { return rotation * Vec3( 1,  0,  0); }
Vec3 Transform::Left()    const { return rotation * Vec3(-1,  0,  0); }
Vec3 Transform::Up()      const { return rotation * Vec3( 0,  1,  0); }
Vec3 Transform::Down()    const { return rotation * Vec3( 0, -1,  0); }
Vec3 Transform::Forward() const { return rotation * Vec3( 0,  0, -1); }
Vec3 Transform::Back()    const { return rotation * Vec3( 0,  0,  1); }

void Transform::Rotate(const Vec3& axis, Float radians) {
  rotation = glm::quat_cast(glm::rotate(Mat4(one), radians, axis)) * rotation;
}

void Transform::RotateAround(const Vec3& axis, const Vec3& pivot, Float radians) {
  position = (glm::quat_cast(glm::rotate(Mat4(one), radians, axis)) * (position - pivot)) + pivot;
}

Vec3 Transform::TransformPoint(const Vec3& point) const {
  Mat4 model = ModelMat();
  return TransformPoint(point, model);
}

/*static*/ Vec3 Transform::TransformPoint(const Vec3& point, const Mat4& model) {
  return { model[0][0] * point.x + model[1][0] * point.y + model[2][0] * point.z + model[3][0],
           model[0][1] * point.x + model[1][1] * point.y + model[2][1] * point.z + model[3][1],
           model[0][2] * point.x + model[1][2] * point.y + model[2][2] * point.z + model[3][2] };
}
