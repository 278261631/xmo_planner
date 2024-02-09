from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='aladin'),
    path('draw_sun', views.draw_sun, name='draw_sun'),
    path('draw_all_sky_by_jgg', views.draw_all_sky_by_jgg, name='draw_all_sky_by_jgg'),
    path('draw_load_plan', views.draw_load_plan, name='draw_load_plan'),
    path('draw_load_jgg_plan', views.draw_load_jgg_plan, name='draw_load_jgg_plan'),
    path('clear_session_data', views.clear_session_data, name='clear_session_data'),
    path('generate_jgg_from_session', views.generate_jgg_from_session, name='generate_jgg_from_session'),
    path('draw_plan', views.draw_plan, name='draw_plan'),
]
