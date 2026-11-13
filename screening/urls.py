from django.urls import path

from . import views

urlpatterns = [
    path("", views.upload, name="upload"),
    path("analysis/<int:pk>/", views.result, name="result"),
    path("analysis/<int:pk>/csv/", views.result_csv, name="result_csv"),
    path("history/", views.history, name="history"),
]
