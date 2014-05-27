from django.conf.urls import patterns
from django.shortcuts import redirect


urlpatterns = patterns('',
    (r"^favicon.ico", "page_not_found"),
    (r"^(?P<path>.+/)$", "viewer.views.index"),
    (r"^(?P<path>.+/).*\.fasta$", "viewer.views.fasta"),
    (r"^$", lambda(x) : redirect("/orthologs/") ),
)
