#ifndef PPCA_SRC_HPP
#define PPCA_SRC_HPP

#include "math.h"
#include <vector>
#include <algorithm>

class Controller {

public:
    Controller(const Vec &_pos_tar, double _v_max, double _r, int _id, Monitor *_monitor) {
        pos_tar = _pos_tar;
        v_max = _v_max;
        r = _r;
        id = _id;
        monitor = _monitor;
    }

    void set_pos_cur(const Vec &_pos_cur) {
        pos_cur = _pos_cur;
    }

    void set_v_cur(const Vec &_v_cur) {
        v_cur = _v_cur;
    }

private:
    int id;
    Vec pos_tar;
    Vec pos_cur;
    Vec v_cur;
    double v_max, r;
    Monitor *monitor;

    bool will_collide_with(const Vec &v_candidate, int other_id) const {
        if (other_id == id) return false;
        Vec pos_j = monitor->get_pos_cur(other_id);
        Vec v_j = monitor->get_v_cur(other_id);
        double r_j = monitor->get_r(other_id);

        Vec delta_pos = pos_cur - pos_j;
        Vec delta_v = v_candidate - v_j;
        double delta_v_norm = delta_v.norm();

        double p = delta_pos.dot(delta_v);
        if (p >= 0) {
            // Not approaching according to the simulator's check
            return false;
        }

        double min_dis_sqr;
        double delta_r = r + r_j;

        if (delta_v_norm < 1e-12) {
            // Relative velocity nearly zero; minimal distance is current
            min_dis_sqr = delta_pos.norm_sqr();
        } else {
            double project = p / (-delta_v_norm);
            if (project < delta_v_norm * TIME_INTERVAL) {
                min_dis_sqr = delta_pos.norm_sqr() - project * project;
            } else {
                min_dis_sqr = (delta_pos + delta_v * TIME_INTERVAL).norm_sqr();
            }
        }
        return min_dis_sqr <= delta_r * delta_r - EPSILON;
    }

    bool is_safe_velocity(const Vec &v_candidate) const {
        int n = monitor->get_robot_number();
        for (int j = 0; j < n; ++j) {
            if (j == id) continue;
            if (will_collide_with(v_candidate, j)) return false;
        }
        return true;
    }

public:

    Vec get_v_next() {
        // If already at target (or extremely close), stop
        Vec to_tar = pos_tar - pos_cur;
        double dist = to_tar.norm();
        if (dist <= EPSILON) {
            return Vec();
        }

        // Desired velocity towards target, limited to v_max and avoiding overshoot
        Vec v_desired;
        double max_step = v_max * TIME_INTERVAL;
        if (dist <= max_step) {
            // Move exactly to the target within this interval
            v_desired = to_tar / TIME_INTERVAL;
        } else {
            v_desired = to_tar.normalize() * v_max;
        }

        // If last step had speeding, be extra conservative (cap to 90%)
        if (monitor->get_speeding(id)) {
            v_desired = v_desired * 0.9;
        }

        // Try to find a safe scaling factor alpha in [0,1]
        // such that no predicted collisions occur (assuming others keep current v)
        double lo = 0.0, hi = 1.0;
        Vec best = Vec();
        // If zero velocity is unsafe (shouldn't be), return zero early
        if (!is_safe_velocity(best)) {
            return Vec();
        }

        // Quick accept if full speed looks safe
        if (is_safe_velocity(v_desired)) {
            return v_desired;
        }

        // Binary search for maximum safe alpha
        for (int it = 0; it < 40; ++it) {
            double mid = (lo + hi) * 0.5;
            Vec v_try = v_desired * mid;
            if (is_safe_velocity(v_try)) {
                best = v_try;
                lo = mid;
            } else {
                hi = mid;
            }
        }

        return best;
    }
};


#endif //PPCA_SRC_HPP

