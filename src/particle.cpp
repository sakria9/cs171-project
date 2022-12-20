#include "particle.h"
#include "object.h"
#include "scene.h"
#include "shader.h"
#include "transform.h"
#include <memory>
#include <map>
#include <vector>
#include <random>
#define INTERACTION_RADIUS 0.05f
#define TEMPATURE 300.0f
#define CORRECTION 1e-28
//密度，光滑核
static Vec3 rand_acc(){
return Vec3((float)rand()/RAND_MAX-0.5f,(float)rand()/RAND_MAX-0.5f,(float)rand()/RAND_MAX-0.5f);
}
static float module(Vec3 vector){
  float ret=(float)sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
  if(ret<1e-4)
  {return 0.0001f;}
  else{
    return ret;
  }
}
static Vec3 normalize(Vec3 vector){
  float mod=module(vector);
  return Vec3(vector.x/mod,vector.y/mod,vector.z/mod);
}
// static float Wpoly6(float r, float h) {
//   if (r > h) return 0;
//   float x = (h * h - r * r) ;
//   return 315.0f / (64.0f * M_PI * pow(h, 9)) * x * x * x;
// }
// //水的刚度系数是2200MPa，这里取2000
// static float WCSPH(float rho, float rho0) {
//   return 0.1f * (pow(rho / rho0, 7) - 1);
// }
// static Vec3 WspikyGrad(float r, Vec3 diff, float h) {
//   if (r > h) return Vec3(0, 0, 0);
//   float x = (h - r) ;
//   return ((float)(45.0f / (M_PI * pow(h, 6)))) *x*x* diff;
// }
// static Vec3 pressureAcc(float rho, float rho0, Vec3 diff, float h) {
//   float coeff=(float)(particle_density
//   //*pow(particle_radius,3)
//   *(WCSPH(rho, rho0)/pow(rho,2)+WCSPH(rho0,rho)/pow(rho0,2)));
//   if (coeff>1e8){
//     //printf("coeff=%f",coeff);
//   }
//   return coeff*
//   WspikyGrad(module(diff),normalize(diff),h);
// }

//重写:
static float Rho(float r, float h) {
  if (r > h) return 0;
  float x = (h * h - r * r) ;
  //printf("%f ",315.0f / (64.0f * M_PI * pow(h, 9)) * x * x * x);
  //from 100 to 1000 value
  float ret=315.0f / (64.0f * M_PI * pow(h, 9)) * x * x * x;
  if(ret<1){
    return 1.0f;
  }
  return ret;
}
static float pressure(float rho, float rho0,float tempature) {
  return tempature*(rho-rho0);
}
static Vec3 pressureAcc(float rhoi, float rhoj,float rho0,float mass,float tempature, Vec3 diff, float h) {
  float coeff=((float)(mass*45/(M_PI*pow(h,6)) 
    *((pressure(rhoi,rho0,tempature)+pressure(rhoj,rho0,tempature))/(rhoi*rhoj)) 
    *((h-module(diff)))
    *((h-module(diff)))));
    if(coeff>100){
    printf("%f ",coeff);}
  Vec3 dir=(normalize(diff));
  if(dir.x!=dir.x||dir.y!=dir.y||dir.z!=dir.z){
    //printf("dir is nan");
    //printf("\ndiff:%f %f %f",diff.x,diff.y,diff.z);
  }
  return coeff * dir;
}
void ParticleSystem::simulate() {
  // TODO: SPH solver !!!!!
  //Note particle_radius=0.01f
  simulateTime++;
  int debugCnt=0;
  std::map<std::vector<int>,std::vector<Particle*>> ptr_map;
  for (auto& particle: particles) {
    std::vector<int> index;
    index.push_back((int)(particle.x[0]/INTERACTION_RADIUS));
    index.push_back((int)(particle.x[1]/INTERACTION_RADIUS));
    index.push_back((int)(particle.x[2]/INTERACTION_RADIUS));
    if(ptr_map.find(index)==ptr_map.end()){
      ptr_map.insert(std::pair<std::vector<int>,std::vector<Particle*>>(index,std::vector<Particle*>()));
      ptr_map[index].push_back(&particle);
      debugCnt++;
    }
    else{
      ptr_map[index].push_back(&particle);
    }
  }
  printf("%d ,", ptr_map.size());
  if(!initialized){
      for (auto& particle: particles) {
    std::vector<int> index;
    index.push_back((int)(particle.x[0]/INTERACTION_RADIUS));
    index.push_back((int)(particle.x[1]/INTERACTION_RADIUS));
    index.push_back((int)(particle.x[2]/INTERACTION_RADIUS));
    Vec3 pressure_acc=Vec3(0,0,0);
    float new_rho=0;
    int x,y,z;
    if((particle.x[0]/INTERACTION_RADIUS)-index[0]>0.5f){
      x=1;
    }
    else{
      x=0;
    }
    if((particle.x[1]/INTERACTION_RADIUS)-index[1]>0.5f){
      y=1;
    }
    else{
      y=0;
    }
    if((particle.x[2]/INTERACTION_RADIUS)-index[2]>0.5f){
      z=1;
    }
    else{
      z=0;
    }

    for(int i=x-1;i<=x;i++){
      for(int j=y-1;j<=y;j++){
        for(int k=z-1;k<=z;k++){
          //std::vector<int> index;
          if(ptr_map.find(index)!=ptr_map.end()){
            for(auto& ptr:ptr_map[std::vector<int>{index[0]+i,index[1]+j,index[2]+k}]){
              if(ptr!=&particle){
                new_rho+=Rho(module(particle.x-ptr->x),INTERACTION_RADIUS);
              }
            }
          }
        }
      }
    }
  particle.rho=new_rho;
  particle.newRho=new_rho;
  }
  initialized=true;
  return;
    }
  for (auto& particle: particles) {
    std::vector<int> index;
    index.push_back((int)(particle.x[0]/INTERACTION_RADIUS));
    index.push_back((int)(particle.x[1]/INTERACTION_RADIUS));
    index.push_back((int)(particle.x[2]/INTERACTION_RADIUS));
    Vec3 pressure_acc=Vec3(0,0,0);
    float new_rho=0;
    int x,y,z;
    if((particle.x[0]/INTERACTION_RADIUS)-index[0]>0.5f){
      x=1;
    }
    else{
      x=0;
    }
    if((particle.x[1]/INTERACTION_RADIUS)-index[1]>0.5f){
      y=1;
    }
    else{
      y=0;
    }
    if((particle.x[2]/INTERACTION_RADIUS)-index[2]>0.5f){
      z=1;
    }
    else{
      z=0;
    }

    for(int i=x-1;i<=x;i++){
      for(int j=y-1;j<=y;j++){
        for(int k=z-1;k<=z;k++){
          //std::vector<int> index;
          if(ptr_map.find(index)!=ptr_map.end()){
            for(auto& ptr:ptr_map[std::vector<int>{index[0]+i,index[1]+j,index[2]+k}]){
              if(ptr!=&particle){
                //debugCnt++;
                // Vec3 diff=particle.x-ptr->x;
                // Float dist=diff.norm();
                // if(dist<INTERACTION_RADIUS){
                //   Vec3 force=diff/dist*(INTERACTION_RADIUS-dist);
                //   particle.a+=force;
                //   ptr->a-=force;
                // }
                //#Example: 斥力
                //new_rho+=particle_density*Wpoly6(module(particle.x-ptr->x),INTERACTION_RADIUS);
                pressure_acc+=pressureAcc(particle.rho,ptr->rho,1.0f,particle.mass,TEMPATURE,particle.x-ptr->x,INTERACTION_RADIUS);
                // if(pressure_acc.x!=pressure_acc.x||pressure_acc.y!=pressure_acc.y||pressure_acc.z!=pressure_acc.z){
                //   printf("%f %f %f ;",pressure_acc.x,pressure_acc.y,pressure_acc.z);
                // }
                if(pressure_acc.x!=pressure_acc.x||pressure_acc.y!=pressure_acc.y||pressure_acc.z!=pressure_acc.z){
                  printf("%f %f %f ;",pressure_acc.x,pressure_acc.y,pressure_acc.z);
                  printf("SimulateTime:%d" ,simulateTime);
                  assert(false);
                }
              }
            }
          }
        }
      }
    }
    //printf("%f %f %f",pressure_acc[0],pressure_acc[1],pressure_acc[2]);
  particle.a=-9.8f*Vec3(0,1,0)+((float)CORRECTION)*pressure_acc;//+(TEMPATURE-273)/100.0f*rand_acc();
  //printf("%f,%f,%f;",pressure_acc.x,pressure_acc.y,pressure_acc.z);
  particle.newRho=new_rho;
  }
  
  for (auto& particle: particles){
    particle.v+=particle.a*fixed_delta_time;
    particle.x+=particle.v*fixed_delta_time;
    particle.rho=particle.newRho;
    clampInBoundaries(particle);
  }
  
  //use isinboundaries to check if the particle is in the box and move it.
  
}

ParticleSystem::ParticleSystem() {}
ParticleSystem::ParticleSystem(int n) {
  boundaries.push_back(Vec3(xmin, ymin, zmin));
  boundaries.push_back(Vec3(xmax, ymin, zmin));
  boundaries.push_back(Vec3(xmin, ymax, zmin));
  boundaries.push_back(Vec3(xmax, ymax, zmin));
  boundaries.push_back(Vec3(xmin, ymin, zmax));
  boundaries.push_back(Vec3(xmax, ymin, zmax));
  boundaries.push_back(Vec3(xmin, ymax, zmax));
  boundaries.push_back(Vec3(xmax, ymax, zmax));

  Vec3 center(0, 1, 1);
  for (int x = -n; x <= n; x++) {
    for (int y = -n; y <= n; y++) {
      for (int z = -n; z <= n; z++) {
        Vec3 pos = center + Vec3(x, y, z) * particle_radius * 2.0f;
        if (isInBoundaries(pos)) {
          particles.push_back(Particle{pos, Vec3(0, 0, 0), Vec3(0, 0, 0)});
        }
      }
    }
  }
}

bool ParticleSystem::isInBoundaries(const Vec3 &x) {
  return xmin <= x.x && x.x <= xmax && ymin <= x.y && x.y <= ymax &&
         zmin <= x.z && x.z <= zmax;
}
void ParticleSystem::clampInBoundaries(Particle &particle) {
  if(particle.x.x<xmin) {particle.x.x=xmin;particle.v.x=0;}
  else if(particle.x.x>xmax) {particle.x.x=xmax;particle.v.x=0;}
  if(particle.x.y<ymin) {particle.x.y=ymin;particle.v.y=0;}
  else if(particle.x.y>ymax) {particle.x.y=ymax;particle.v.y=0;}
  if(particle.x.z<zmin) {particle.x.z=zmin;particle.v.z=0;}
  else if(particle.x.z>zmax) {particle.x.z=zmax;particle.v.z=0;}
}

std::vector<std::shared_ptr<Object>> ParticleSystem::boundryIndicators() {
  const Float scale = 0.05;
  std::vector<std::shared_ptr<Object>> ret;
  for (auto &x : boundaries) {
    auto obj = std::make_shared<Object>(
        mesh_sphere, Shader::shader_phong,
        Transform(x, Quat(1, 0, 0, 0), Vec3(scale, scale, scale)));
    obj->color = {one, 0, 0};
    ret.emplace_back(obj);
  }
  return ret;
}

void ParticleSystem::simulate(unsigned int num_steps) {
  for (unsigned i = 0; i < num_steps; ++i)
    simulate();
}

void ParticleSystem::fixedUpdate() {
  simulate(simulation_steps_per_fixed_update_time);
}

void ParticleSystem::renderParticle(const Scene &scene) {
  const Float particle_scale = 0.8;
  auto shader = Shader::shader_phong;
  const Vec3 color(0, 0, one);
  auto transform = Transform(Vec3(0, 0, 0), Quat(1, 0, 0, 0),
                             Vec3(particle_radius * 2 * particle_scale,
                                  particle_radius * 2 * particle_scale,
                                  particle_radius * 2 * particle_scale));
  for (auto &particle : particles) {
    transform.position = particle.x;
    shader->Set("model", transform.ModelMat());
    shader->Set("view", scene.camera.LookAtMat());
    shader->Set("projection", scene.camera.PerspectiveMat());
    shader->Set("object_color", color);
    shader->Set("light_position", scene.light_position);
    shader->Set("light_color", scene.light_color);
    shader->Set("camera_position", scene.camera.transform.position);
    shader->Set("is_bidirectional", false);
    shader->Use();
    mesh_sphere->DrawTriangles();
  }
}
