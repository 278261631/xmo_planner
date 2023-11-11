import astropy.units as u
import numpy as np
from astropy.coordinates import AltAz, GCRS, SkyCoord


def rotate(np_array_a, or_theta_rad, or_axis):
    or_axis = or_axis / np.linalg.norm(or_axis)  # 确保是单位向量
    # 构建旋转矩阵
    c_rotate = np.cos(or_theta_rad)
    s_rotate = np.sin(or_theta_rad)
    t_rotate = 1 - c_rotate
    x_rotate, y_rotate, z_rotate = or_axis
    r_rotate = np.array([
        [t_rotate * x_rotate * x_rotate + c_rotate, t_rotate * x_rotate * y_rotate - z_rotate * s_rotate,
         t_rotate * x_rotate * z_rotate + y_rotate * s_rotate],
        [t_rotate * x_rotate * y_rotate + z_rotate * s_rotate, t_rotate * y_rotate * y_rotate + c_rotate,
         t_rotate * y_rotate * z_rotate - x_rotate * s_rotate],
        [t_rotate * x_rotate * z_rotate - y_rotate * s_rotate, t_rotate * y_rotate * z_rotate + x_rotate * s_rotate,
         t_rotate * z_rotate * z_rotate + c_rotate]
    ])
    # 将向量 np_array_a 与旋转矩阵相乘
    a_rotated = np.dot(r_rotate, np_array_a)
    print("原始向量 np_array_a:", np_array_a)
    print("旋转后的向量 a_rotated:", a_rotated)
    return a_rotated


def get_intersection_point(point_pa, point_pb, point_pc, point_pd):
    # 向量BC和BD
    vector_bc = point_pc - point_pb
    vector_bd = point_pd - point_pb
    # # 计算平面BCD的法向量
    normal_vector_bcd = np.cross(vector_bc, vector_bd)
    # 求平面BCD的方程，ax + by + cz + d_val = 0
    ax, by, cz = normal_vector_bcd
    d_val = -np.dot(normal_vector_bcd, point_pb)
    # 求A到平面BCD的垂线与平面BCD的交点坐标
    distance_to_plane = (ax * point_pa[0] + by * point_pa[1] + cz * point_pa[2] + d_val) / (ax ** 2 + by ** 2 + cz ** 2)
    intersection_point = point_pa - distance_to_plane * normal_vector_bcd
    # 输出交点坐标
    print("交点坐标:", intersection_point)
    return intersection_point


def get_ground_plane_intersection(point_x, obs_time=None, obs_location=None):
    az_north = 0.0 * u.deg
    az_east = 90 * u.deg
    alt_temp = 0.0 * u.deg
    # coord_north = SkyCoord(alt=alt_temp, az=az_north, frame='altaz', obstime=obs_time, location=obs_location)
    # coord_east = SkyCoord(alt=alt_temp, az=az_east, frame='altaz', obstime=obs_time, location=obs_location)
    # coord_north = coord_north.transform_to(AltAz(obs_time, obs_location)).transform_to('gcrs')
    # coord_east = coord_east.transform_to(AltAz(obs_time, obs_location)).transform_to('gcrs')

    # 创建AltAz坐标
    coord_north = AltAz(az=az_north, alt=alt_temp, obstime=obs_time, location=obs_location)
    coord_east = AltAz(az=az_east, alt=alt_temp, obstime=obs_time, location=obs_location)

    coord_north = coord_north.transform_to(GCRS)
    coord_east = coord_east.transform_to(GCRS)

    plan_center = np.array([0, 0, 0])
    north_cart_array = np.array(
        [coord_north.cartesian.x.value, coord_north.cartesian.y.value, coord_north.cartesian.z.value])
    east_cart_array = np.array(
        [coord_east.cartesian.x.value, coord_east.cartesian.y.value, coord_east.cartesian.z.value])
    return get_intersection_point(point_x, north_cart_array, east_cart_array, plan_center)


def get_star_plane_intersection(point_star_x, obs_time=None, obs_location=None):
    point_ground_x = get_ground_plane_intersection(point_star_x, obs_time, obs_location)
    ax, by, cz = point_star_x
    plan_center = np.array([0, 0, 0])
    d_val = -np.dot(point_star_x, plan_center)
    # 求A到平面BCD的垂线与平面BCD的交点坐标
    distance_to_plane = (ax * point_ground_x[0] + by * point_ground_x[1] + cz * point_ground_x[2] + d_val) / (
                ax ** 2 + by ** 2 + cz ** 2)
    intersection_point = point_ground_x - distance_to_plane * point_star_x
    return intersection_point


def get_star_plane_intersection_by_z_axis(point_star_x):
    z_axis_point = np.array([0, 0, 1])
    ax, by, cz = point_star_x
    plan_center = np.array([0, 0, 0])
    d_val = -np.dot(point_star_x, plan_center)
    # 求A到平面BCD的垂线与平面BCD的交点坐标
    distance_to_plane = (ax * z_axis_point[0] + by * z_axis_point[1] + cz * z_axis_point[2] + d_val) / (
                ax ** 2 + by ** 2 + cz ** 2)
    intersection_point = z_axis_point - distance_to_plane * point_star_x
    return intersection_point


def get_l_r_t_b_axis(point_star_x, obs_time=None, obs_location=None):
    intersection_point = get_star_plane_intersection(point_star_x, obs_time, obs_location)
    # left_r_axis = np.linalg.norm(intersection_point)
    left_r_axis = intersection_point / np.linalg.norm(intersection_point)
    right_r_axis = np.negative(left_r_axis)
    top_r_axis = rotate(left_r_axis, np.pi / 2, point_star_x)
    # bottom_axis = rotate(left_axis, np.pi/2, point_star_x)
    bottom_r_axis = np.negative(top_r_axis)
    return left_r_axis, right_r_axis, top_r_axis, bottom_r_axis


def get_l_r_t_b_axis_by_az(point_star_x):
    intersection_point = get_star_plane_intersection_by_z_axis(point_star_x)
    left_r_axis = intersection_point / np.linalg.norm(intersection_point)
    right_r_axis = np.negative(left_r_axis)
    top_r_axis = rotate(left_r_axis, np.pi / 2, point_star_x)
    # bottom_axis = rotate(left_axis, np.pi/2, point_star_x)
    bottom_r_axis = np.negative(top_r_axis)
    print('l:%s  (%s,%s,%s)' % (left_r_axis, left_r_axis[0], left_r_axis[1], left_r_axis[2]))
    print('r:%s  (%s,%s,%s)' % (right_r_axis, right_r_axis[0], right_r_axis[1], right_r_axis[2]))
    print('t:%s  (%s,%s,%s)' % (top_r_axis, top_r_axis[0], top_r_axis[1], top_r_axis[2]))
    print('b:%s  (%s,%s,%s)' % (bottom_r_axis, bottom_r_axis[0], bottom_r_axis[1], bottom_r_axis[2]))
    return left_r_axis, right_r_axis, top_r_axis, bottom_r_axis


def get_left(c_ra_deg, c_dec_deg, theta_left_deg, obs_time=None, obs_location=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    # left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis(center_cart_array, obs_time, obs_location)
    left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis_by_az(center_cart_array)
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_left_deg), left_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_right(c_ra_deg, c_dec_deg, theta_right, obs_time=None, obs_location=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    # left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis(center_cart_array, obs_time, obs_location)
    left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis_by_az(center_cart_array)
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_right), right_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_top(c_ra_deg, c_dec_deg, theta_top, obs_time=None, obs_location=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    # left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis(center_cart_array, obs_time, obs_location)
    left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis_by_az(center_cart_array)
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_top), top_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_bottom(c_ra_deg, c_dec_deg, theta_bottom, obs_time=None, obs_location=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    # left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis(center_cart_array, obs_time, obs_location)
    left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis_by_az(center_cart_array)
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_bottom), bottom_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_rotate_fix_axis(c_ra_deg, c_dec_deg):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    left_axis_c, right_axis_c, top_axis_c, bottom_axis_c = get_l_r_t_b_axis_by_az(center_cart_array)

    return left_axis_c, right_axis_c, top_axis_c, bottom_axis_c


def get_left_fix_axis(c_ra_deg, c_dec_deg, theta_left_deg, left_axis_c=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_left_deg), left_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_right_fix_axis(c_ra_deg, c_dec_deg, theta_right, right_axis_c=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_right), right_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_top_fix_axis(c_ra_deg, c_dec_deg, theta_top, top_axis_c=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_top), top_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')


def get_bottom_fix_axis(c_ra_deg, c_dec_deg, theta_bottom, bottom_axis_c=None):
    center_coord = SkyCoord(ra=c_ra_deg, dec=c_dec_deg, unit='deg')
    center_cart = center_coord.cartesian
    print('cart  [%s]' % center_cart)
    center_cart_array = np.array([center_cart.x.value, center_cart.y.value, center_cart.z.value])
    center_cart_l = rotate(center_cart_array, np.deg2rad(theta_bottom), bottom_axis_c)
    cen_cart_coord_l = SkyCoord(x=center_cart_l[0] * u.pc, y=center_cart_l[1] * u.pc, z=center_cart_l[2] * u.pc,
                                representation_type='cartesian')
    return cen_cart_coord_l.transform_to('gcrs')
