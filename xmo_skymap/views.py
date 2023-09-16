import pytz

from astropy.coordinates import EarthLocation, AltAz, get_sun, get_moon, Angle, get_body
from astropy.time import Time
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from datetime import datetime
from django.core.cache import cache
from django.template import loader
import astropy.units as u


def index(request):

    # return HttpResponse("skymap")
    now_date = datetime.today().strftime('%Y%m%d')
    req_date = request.GET.get('req_date', now_date)
    print(' -> req [%s]' % req_date)

    latest_tele_system_list = []

    print(cache.get(req_date))
    sky_process_map = cache.get(req_date)
    template = loader.get_template('skymap/aladin.html')
    context = {
        'latest_tele_system_list': latest_tele_system_list,
        # 'day_in_cache': num,
        'process_map': sky_process_map,
        'req_date': req_date,

    }
    return HttpResponse(template.render(context, request))


def draw_sun(request):
    now_date = datetime.today().strftime('%Y%m%d')
    url = request.GET.get('url_for_test', now_date)
    ex_message = ''

    # 设置观测地点（以纽约为例）
    location = EarthLocation(lat='39.7128', lon='116.0060')  # 纬度和经度

    # 获取当前时间
    current_time = Time.now()
    current_time_with_tz = current_time.to_datetime(pytz.timezone('Asia/Shanghai'))
    # 将观测地点和时间结合起来
    altaz = AltAz(location=location, obstime=current_time)

    # 计算太阳的赤经和赤纬
    sun = get_body('sun', current_time, location)
    print(sun)
    sun_radec = sun.transform_to(AltAz(location=location, obstime=current_time)).transform_to('gcrs')

    # 计算月亮的赤经和赤纬
    moon = get_body('moon', current_time, location)
    moon_radec = moon.transform_to(AltAz(location=location, obstime=current_time)).transform_to('gcrs')
    print(moon)
    sun_ra_hms = Angle(sun_radec.ra.deg, unit=u.deg)
    moon_ra_hms = Angle(moon_radec.ra.deg, unit=u.deg)

    # ra_hms_str = ra_hms.to_string(unit=u.hour, sep=':', precision=2)
    # 打印太阳和月亮的赤经和赤纬
    print(f"太阳的赤经（Right Ascension）：{sun_radec.ra.deg}  {sun_ra_hms.hms}  {sun_ra_hms.dms}")
    print(f"太阳的赤纬（Declination）：{sun_radec.dec.deg}")
    print()
    print(f"月亮的赤经（Right Ascension）：{moon_radec.ra.deg}   {moon_ra_hms.hms}   {moon_ra_hms.dms}")
    print(f"月亮的赤纬（Declination）：{moon_radec.dec.deg}")
    print(f"time {current_time.iso}")
    print(f"当前时间：{current_time_with_tz.strftime('%Y-%m-%d %H:%M:%S %Z')}")

    return JsonResponse({'sun_ra': str(sun_radec.ra.deg), 'sun_dec': str(sun_radec.dec.deg),
                         'moon_ra': str(moon_radec.ra.deg), 'moon_dec': str(moon_radec.dec.deg),
                         'ex_message': ex_message})

