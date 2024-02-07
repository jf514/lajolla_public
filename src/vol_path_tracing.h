#pragma once

#include <sstream>

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources

Vector3 vecExp(const Vector3& vec){
    return Vector3( std::exp(vec.x), std::exp(vec.y), std::exp(vec.z));
}

extern bool gdb;

static std::mutex printMut;
void printdb(const std::string& s1, const double& s2){

    if(!gdb)
        return;

    std::unique_lock<std::mutex> lock(printMut);
    std::cout << s1 << " " << s2 << "\n";
    lock.unlock();
}

void printdb(const std::string& s1, const Vector3& s2){

    if(!gdb)
        return;

    std::unique_lock<std::mutex> lock(printMut);
    std::cout << s1 << ": " << s2.x << ", " << s2.y << ", " << s2.z << "\n";
    lock.unlock();
}

Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{0,0};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if(vertex_){
        PathVertex vertex = *vertex_;
        Spectrum sigma_a = get_sigma_a(scene.media[vertex.exterior_medium_id], vertex.position);
        Vector3 Le{0,0,0};
        if (is_light(scene.shapes[vertex.shape_id])) {
            Spectrum Le = emission(vertex, -ray.dir, scene);
            Real t = dot(vertex.position - ray.org, ray.dir);

            Real transmittance = std::exp(-sigma_a.x * t);
            printdb("out ", transmittance*Le);
            return Le*transmittance;
        }
    }
    
    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {

  //std::unique_lock<std::mutex> lock(printMut);

    //using ss = std::stringstream;

    //Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{0,0};
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    //printdb("num lights", scene.lights.size());

    Real t_hit = infinity<Real>(); // Maybe set to something large but finite
    if(vertex_.has_value()){
        //printdb("here vtx", 0);
        t_hit = dot(ray.dir, vertex_->position - ray.org);
        //printdb("t_hit", t_hit);
    }

    int medium_id = 0;
    Spectrum sigma_s = get_sigma_s(scene.media[medium_id], Vector3{1,1,1});
    printdb("sigma_s", sigma_s);
    Spectrum sigma_a = get_sigma_a(scene.media[medium_id], Vector3{1,1,1});
    Spectrum sigma_t = sigma_a + sigma_s;

    Real u = next_pcg32_real<Real>(rng);
    Real t = - std::log(1-u)/sigma_t[0];

    //printdb("t-thit", t - t_hit);
    //printdb("t ", t);
    //printdb("t_hit", t_hit);

    if(t < t_hit){
        printdb("t<t_hit",0);
        Real trans_pdf = std::exp(-sigma_t[0] * t) * sigma_t[0];
        Real transmittance = std::exp(-sigma_t[0] * t);
        Vector3 p = ray.org + t * ray.dir;
        PhaseFunction phase_func = get_phase_function(scene.media[medium_id]);

        // Next event...
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        //printdb("light id", light_id);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);

        Real dist_to_p = distance(point_on_light.position, p);
        Vector3 dir_light = normalize(point_on_light.position - p);

        Ray shadow_ray{p, dir_light,
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                            distance(point_on_light.position, p)};
        if (!occluded(scene, shadow_ray)){
            printdb("!occluded", 0);

            Spectrum integrand{1,1,1};
            Real pdf = 1;

            // Phase
            integrand *= eval(phase_func, dir_light, ray.dir);
            //pdf *= pdf_sample_phase(phase_func, dir_light, ray.dir);

            // L_e
            Spectrum Le = emission(light, -dir_light, Real(0), point_on_light, scene);
            Real L_pdf = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, p, scene);

            integrand *= Le;
            pdf *= L_pdf;

            Real fac = fabs(dot(dir_light, point_on_light.normal)) * 
                std::exp(-sigma_t[0] * dist_to_p) / (dist_to_p * dist_to_p);

            Vector3 out = (transmittance / trans_pdf) * sigma_s * (integrand / pdf ) * fac;
            //Vector3 out = transmittance * fac * L;
            printdb("out", out);
            return out;
        } 
    } else {
        //printdb("here", 1);

        if(!vertex_.has_value()){
            //printdb("exit", 2);
            return make_zero_spectrum();
        }

        PathVertex vertex = *vertex_;
        printdb("vertex.pos", vertex.position);
        Real trans_pdf = std::exp(-sigma_t[0] * t_hit) * sigma_t[0];
        printdb("t_hit 0", t_hit);
        Real transmittance = std::exp(-sigma_t[0] * t_hit);
 

        Vector3 Le{0,0,0};
        if (is_light(scene.shapes[vertex.shape_id])) {
            Spectrum Le = emission(vertex, -ray.dir, scene);

            // printdb("Le ", Le);
            // printdb("t_hit", t_hit);
            // printdb("transmittance", transmittance);
            // printdb("out", Le*transmittance);
            // printdb("sigma_a", sigma_a);
            // printdb("sigma_t", sigma_t);
            // printdb("sigma_s", sigma_s);
            return Le*transmittance / trans_pdf;
        }        
    } 

    //lock.unlock();
    
    //return make_zero_spectrum();
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
