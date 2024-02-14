//pragma once

#include <sstream>

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources

Vector3 vecExp(const Vector3& vec){
    return Vector3( std::exp(vec.x), std::exp(vec.y), std::exp(vec.z));
}

extern bool gdb;

static std::mutex printMut;
void printdb(const std::string& s1){
   if(!gdb)
        return;

    std::unique_lock<std::mutex> lock(printMut);
    std::cout << s1 << "\n";
    lock.unlock();
}

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

void printdb(const Ray& ray){
    if(!gdb)
        return;

    printdb("ray.org", ray.org);
    printdb("ray.dir", ray.dir);
}

void printdb(const std::string& str, const Ray& ray){
    if(!gdb)
        return;

    std::cout << str << "\n";
    printdb(ray);
}

Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    ray.dir = Vector3(0,0,1);
    //printdb("ray.org", ray.org);
    RayDifferential ray_diff = RayDifferential{0,0};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    printdb("hit loc", vertex_->position);
    if(vertex_){
        PathVertex vertex = *vertex_;
        Spectrum sigma_a = get_sigma_a(scene.media[vertex.exterior_medium_id], vertex.position);
        Vector3 Le{0,0,0};
        if (is_light(scene.shapes[vertex.shape_id])) {
            Spectrum Le = emission(vertex, -ray.dir, scene);
            printdb("Le_si", Le);
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
        ray.dir = Vector3{0,0,1};
        PhaseFunction phase_func = get_phase_function(scene.media[medium_id]);

        // Next event...
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        assert(light_id >= 0);
        //printdb("light id", light_id);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);

        Real dist_to_p = distance(point_on_light.position, p);
        Vector3 dir_to_light = normalize(point_on_light.position - p);

        Ray shadow_ray{p, dir_to_light,
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                            distance(point_on_light.position, p)};
        if (!occluded(scene, shadow_ray)){
            printdb("!occluded", 0);

            Spectrum integrand{1,1,1};
            Real pdf = 1;

            // Phase
            integrand *= eval(phase_func, dir_to_light, ray.dir);
            //pdf *= pdf_sample_phase(phase_func, dir_to_light, ray.dir);

            // L_e
            Spectrum Le = emission(light, -dir_to_light, Real(0), point_on_light, scene);
            Real L_pdf = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, p, scene);

            integrand *= Le;
            pdf *= L_pdf;

            Real fac = fabs(dot(dir_to_light, point_on_light.normal)) * 
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
        Real trans_pdf = std::exp(-sigma_t[0] * t_hit);
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
    assert(false);
    //return make_zero_spectrum();
}

int update_medium(const PathVertex& isect, const Ray& ray, int medium_id);

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {

    //printdb("rr", scene.options.rr_depth);
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    //ray.dir = Vector3(0,0,1);
    printdb("Top Ray:");
    printdb(ray);
    RayDifferential ray_diff = RayDifferential{0,0};


    int current_medium_id = scene.camera.medium_id;
    Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], Vector3{1,1,1});
    printdb("sigma_s", sigma_s);
    Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id], Vector3{1,1,1});
    Spectrum sigma_t = sigma_a + sigma_s;
    printdb("sigma_t", sigma_t);

    Spectrum current_path_throughput{1,1,1};
    Spectrum radiance{0,0,0};
    int bounces = 0;
    printdb("\n\n");
    printdb("************************");
    while (true) {
        //printdb("Bounce", bounces);
        bool scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);
        printdb("\nray.org", ray.org);
        //printdb("ray.dir", ray.dir);
        printdb("sigma = ", Vector3(sigma_a[0], sigma_s[0], sigma_t[0]));

        Real t_hit = infinity<Real>(); 
        if(isect_.has_value()){
            //printdb("HIT");
            //printdb("Surf id ", isect_->shape_id);
            //printdb("here vtx", 0);
            t_hit = dot(ray.dir, isect_->position - ray.org);
            //printdb("t_hit", t_hit);
        } 

        // isect might not intersect a surface, but we might be in a volume
        Spectrum transmittance{1,1,1};
        Real trans_pdf = 1.0;
        if (current_medium_id != -1) {
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], Vector3{1,1,1});
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium_id], Vector3{1,1,1});
            Spectrum sigma_t = sigma_a + sigma_s;
            printdb("check t_hit");
            printdb("t_hit", t_hit);
            printdb("hit loc", isect_->position);
            printdb("sigma = ", Vector3(sigma_a[0], sigma_s[0], sigma_t[0]));

            Real u = next_pcg32_real<Real>(rng);
            Real t = - std::log(1-u)/sigma_t[0];

            // compute transmittance and trans_pdf

            // if t < t_hit, set scatter = True
            //t = t_hit;
            if(t < t_hit){
                transmittance *= std::exp(-sigma_t[0] * t);
                trans_pdf = sigma_t[0] * std::exp(-sigma_t[0] * t);
                scatter = true;
                ray.org = ray.org + (t * ray.dir);
            } else {
                transmittance *= std::exp(-sigma_t[0] * t_hit);
                trans_pdf =  std::exp(-sigma_t[0] * t_hit);
                //trans_pdf = 1;
                printdb("trans_pdf", trans_pdf);
                //printdb("sigma_t", sigma_t[0]);
                
                ray.org = ray.org + (t_hit) * ray.dir;
                printdb("ray.org", ray.org);
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);
        printdb("ctp", 0);
        printdb("curr_tp", current_path_throughput);
        printdb("trans", transmittance);
        printdb("trans_pdf", trans_pdf);

        if (!scatter) {
            // reach a surface, include emission
            Vector3 Le_isect{0,0,0};
            if (isect_ && is_light(scene.shapes[isect_->shape_id])) {
                printdb("Hit light");
                Le_isect = emission(*isect_, -ray.dir, scene);
                printdb("Le_isect", Le_isect);
            }
            radiance += current_path_throughput * Le_isect;
            printdb("hit - radiance", radiance);
            printdb("curr_tp", current_path_throughput);
            //printdb("Le isect", Le_isect);
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            // reach maximum bounces
            break;
        }

        if (!scatter && isect_) {
            if (isect_->material_id == -1) {
                // index-matching interface, skip through it
                printdb("pre-sig update");
                printdb("sigma = ", Vector3(sigma_a[0], sigma_s[0], sigma_t[0]));
                current_medium_id = update_medium(*isect_, ray, current_medium_id);
                sigma_s = get_sigma_s(scene.media[current_medium_id], Vector3{1,1,1});
                sigma_a = get_sigma_a(scene.media[current_medium_id], Vector3{1,1,1});
                sigma_t = sigma_a + sigma_s;
                printdb("post-sig update");
                printdb("sigma = ", Vector3(sigma_a[0], sigma_s[0], sigma_t[0]));

                bounces += 1;
 
                continue;
            }
        }

        // sample next direct & update path throughput
        if (scatter) {
            PhaseFunction phase_func = get_phase_function(scene.media[current_medium_id]);
            Vector2 scatter_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Vector3 next_dir = sample_phase_function(phase_func, -ray.dir, scatter_uv).value();
            // update ray.dir
            ray.dir = next_dir;
            printdb("scatter - next_dir", next_dir);
            printdb("scatter !!!!!");
            //ray.dir = Vector3(0,0,1);
            current_path_throughput *= (eval(phase_func, -ray.dir, next_dir) / pdf_sample_phase(phase_func,-ray.dir, next_dir)) * sigma_s;
        } else {
            // Hit a surface -- don’t need to deal with this yet
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            //printdb("rr", scene.options.rr_depth);
            rr_prob = std::min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
                printdb("ctp",1);
                printdb("curr_tp", current_path_throughput);
            }
        }

        bounces += 1;
        printdb("");
    }
    printdb("radiance - final", radiance);
    return radiance;
}

struct Sigma {
    Real a;
    Real s;
    Real t;
};

Sigma get_sigma(const Scene& scene, int medium_id){
    Sigma sigma;
    if(medium_id != -1){
        const Vector3 UNUSED{1,1,1};
        sigma.s = get_sigma_s(scene.media[medium_id], UNUSED)[0];
        sigma.a = get_sigma_a(scene.media[medium_id], UNUSED)[0];
        sigma.t = sigma.a + sigma.s;
    } else {
        sigma.s = sigma.s = sigma.t = 0;
    }
    return sigma;
}

void CheckLight(const std::string& st, PathVertex& v, const Scene& scene){

    int light_id = get_area_light_id(scene.shapes[v.shape_id]);
    if(light_id == -1){
        std::unique_lock<std::mutex> lock(printMut);
        std::cout << "st";
        lock.unlock();
    }
}

// Result for next event estimation
struct NeeResults {
    Spectrum radiance = Spectrum(0,0,0);
    Real pdf = 0;
    int light_id = -1;
};

NeeResults next_event_estimation(Vector3 p, const Scene& scene, int current_medium_id, 
    
    // Sample light
    int bounces, pcg32_state& rng) {
    Vector3 nee_p_cache = p;
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    assert(light_id >= 0);
    const Light &light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 dir_to_light = normalize(p_prime.position - p);

    // Compute transmittance to light. Skip through index-matching shapes.
    Real T_light = 1;
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Real p_trans_dir = 1; // for multiple importance sampling
    while (true) {
        // Cast shadow ray.
        Ray shadow_ray{p, p_prime.position, 
                        get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                            distance(p_prime.position, p)};

        std::optional<PathVertex> isect_ = intersect(scene, shadow_ray);
        Real next_t = distance(p, p_prime.position);
        if (isect_) {
            next_t = distance(p, isect_->position);
        }
        // Account for the transmittance to next_t
        if (shadow_medium_id != -1) {
            Sigma sigma = get_sigma(scene, shadow_medium_id);
            T_light *= exp(-sigma.t * next_t);
            p_trans_dir *= exp(-sigma.t * next_t);
        }
        if (!isect_) {
            // Nothing is blocking, we’re done
            break;
        } else {
            // Something is blocking: is it an opaque surface?
            if (isect_->material_id >= 0) {
                // we’re blocked
                return NeeResults{};
            }
            // otherwise, it’s an index-matching surface and
            // we want to pass through -- this introduces
            // one extra connection vertex
            shadow_bounces += 1;
            if (scene.options.max_depth != -1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                // Reach the max no. of vertices
                return NeeResults{};
            }
            shadow_medium_id = update_medium(*isect_, shadow_ray, shadow_medium_id);
            p = p + next_t * dir_to_light;
        }
        if (T_light > 0) {
            // Compute T_light * G * rho * L & pdf_nee
            
            // G:
            PhaseFunction phase_func = get_phase_function(scene.media[shadow_medium_id]);
            Vector2 scatter_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Vector3 next_dir = sample_phase_function(phase_func, dir_to_light, scatter_uv).value();
            Real dist_to_p_prime = distance(p, p_prime.position);
            Real G = fabs(dot(dir_to_light, p_prime.normal)) /
                (dist_to_p_prime * dist_to_p_prime);

            // rho
            Spectrum rho = eval(phase_func, -shadow_ray.dir, next_dir);
            
            // L
            Spectrum L = emission(light, -dir_to_light, Real(0), p_prime, scene);

            // pdf_nee
            Real pdf_nee = pdf_point_on_light(light, p_prime, nee_p_cache, scene);
            
            Spectrum contrib = T_light * G * rho * L / pdf_nee;

            // Multiple importance sampling: it’s also possible
            // that a phase function sampling + multiple exponential sampling
            // will reach the light source.
            // We also need to multiply with G to convert phase function PDF to area measure.
            PhaseFunction phase_function;
            Real pdf_phase = pdf_sample_phase(phase_function, dir_to_light, p_prime.normal) * G * p_trans_dir;
            // power heuristics
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            
            NeeResults nee_res;
            nee_res.radiance = w*contrib;
            nee_res.pdf = pdf_nee;
            nee_res.light_id = light_id;
            return nee_res;
        }
    }
    return NeeResults();
}

int update_medium(const PathVertex& isect, const Ray& ray, int medium_id) {
    int out_medium_id = medium_id;
    if (isect.interior_medium_id != isect.exterior_medium_id) {
        // At medium transition. Update medium.
        if (dot(ray.dir, isect.geometric_normal) > 0) {
            out_medium_id = isect.exterior_medium_id;
        } else {
            out_medium_id = isect.interior_medium_id;
        }
    }
    return out_medium_id;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {

    // Cast ray                                    
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{0,0};


    int current_medium_id = scene.camera.medium_id;
    Sigma sigma = get_sigma(scene, current_medium_id);

    Spectrum current_path_throughput{1,1,1};
    Spectrum radiance{0,0,0};
    int bounces = 0;

    // Variables for loop invariant
    Real dir_pdf = 1;
    Vector3 nee_p_cache{0,0,0};
    Real multi_trans_pdf = 1;
    bool scatter = false;
    bool never_scattered = true;
    Real pdf_nee = 0; 
    while (true) {
        if(scatter){
            // Debug note - we are making it here!
            never_scattered = false;
        }
        scatter = false;
        std::optional<PathVertex> isect_ = intersect(scene, ray, ray_diff);

        Real t_hit = infinity<Real>(); 
        if(isect_.has_value()){
            t_hit = dot(ray.dir, isect_->position - ray.org);
        } 
        // isect might not intersect a surface, but we might be in a volume
        Spectrum transmittance{1,1,1};
        Real trans_pdf = 1.0;
        if (current_medium_id != -1) {
            Sigma sigma = get_sigma(scene, current_medium_id);

            Real u = next_pcg32_real<Real>(rng);
            Real t = - std::log(1-u)/sigma.t;

            // compute transmittance and trans_pdf
            if(t < t_hit){
                transmittance *= std::exp(-sigma.t * t);
                trans_pdf = sigma.t * std::exp(-sigma.t * t);
                scatter = true;
                ray.org = ray.org + (t * ray.dir);
            }
        } else if(isect_.has_value()) {
            transmittance *= std::exp(-sigma.t * t_hit);
            trans_pdf =  std::exp(-sigma.t * t_hit);
            constexpr Real epsilon = 0.001;
            ray.org = ray.org + t_hit * ray.dir - epsilon * isect_->geometric_normal;
        }
        
        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= (transmittance / trans_pdf);
        if (!scatter) {    
            if(never_scattered){
                // This is the only way we can see the light source, so
                // we don’t need multiple importance sampling.
                Vector3 Le_isect{0,0,0};
                if (isect_ && is_light(scene.shapes[isect_->shape_id])) {
                    Le_isect = emission(*isect_, -ray.dir, scene);
                }
                radiance += current_path_throughput * Le_isect;
            } else if(is_light(scene.shapes[isect_->shape_id])){
                PointAndNormal light_point{isect_->position, isect_->geometric_normal};
                //Need to account for next event estimation
                //Note that pdf_nee needs to account for the path vertex that issued
                //next event estimation potentially many bounces ago.
                //The vertex position is stored in nee_p_cache.
                const int light_id = get_area_light_id(scene.shapes[isect_->shape_id]);
                const Light& light = scene.lights[light_id];
                Real pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene);

                //The PDF for sampling the light source using phase function sampling + transmittance sampling
                //The directional sampling pdf was cached in dir_pdf in solid angle measure.
                //The transmittance sampling pdf was cached in multi_trans_pdf.
                Real G = fabs(dot(-ray.dir, isect_->geometric_normal));
                Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                //current_path_throughput already accounts for transmittance.
                radiance += w * current_path_throughput * emission(*isect_, -ray.dir, scene);
            }
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1) {
            // reach maximum bounces
            break;
        }

        if (!scatter && isect_) {
            if (isect_->material_id == -1) {
                // index-matching interface, skip through it
                current_medium_id = update_medium(*isect_, ray, current_medium_id);
                sigma = get_sigma(scene, current_medium_id);
                bounces += 1;
 
                continue;
            }
        }

        // sample next direct & update path throughput
        if (scatter) {
            PhaseFunction phase_func = get_phase_function(scene.media[current_medium_id]);
            Vector2 scatter_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Vector3 next_dir = sample_phase_function(phase_func, -ray.dir, scatter_uv).value();
            // update ray.dir
            ray.dir = next_dir;
            Real phase_pdf = pdf_sample_phase(phase_func,-ray.dir, next_dir);
            current_path_throughput *= (eval(phase_func, -ray.dir, next_dir) / phase_pdf) * sigma.s;
            dir_pdf *= phase_pdf;
            Real w = (dir_pdf * dir_pdf) / (dir_pdf * dir_pdf + pdf_nee * pdf_nee);
            NeeResults res = next_event_estimation(ray.org, scene, current_medium_id, bounces, rng);
            if(res.light_id != -1){
                radiance += w * res.radiance;
                pdf_nee = res.pdf;
            }
            
        } else {
            // Hit a surface -- don’t need to deal with this yet
            break;
        }

        // Russian roulette
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = std::min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }

        bounces += 1;
    }
    return radiance;
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
