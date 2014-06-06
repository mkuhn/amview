# Django settings for amview project.

import os

DEBUG = not False
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    ('Michael Kuhn', 'michael.kuhn@gmail.com'),
)

MANAGERS = ADMINS

ALLOWED_HOSTS = [ 'djangosrv', '.tu-dresden.de' ]

if not DEBUG:
    CACHES = {
        'default': {
            'BACKEND': 'django.core.cache.backends.filebased.FileBasedCache',
            'LOCATION': '/tmp/django_cache',
            'OPTIONS': {
                'MAX_ENTRIES': 10000
            },
            'VERSION' : 3,
        }
    }


DATABASE_ENGINE = ''           # 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
DATABASE_NAME = ''             # Or path to database file if using sqlite3.
DATABASE_USER = ''             # Not used with sqlite3.
DATABASE_PASSWORD = ''         # Not used with sqlite3.
DATABASE_HOST = ''             # Set to empty string for localhost. Not used with sqlite3.
DATABASE_PORT = ''             # Set to empty string for default. Not used with sqlite3.

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'America/Chicago'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = False

# Absolute path to the directory that holds media.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = ''

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash if there is a path component (optional in other cases).
# Examples: "http://media.lawrence.com", "http://example.com/media/"
MEDIA_URL = ''

# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = '/media/'

# Make this unique, and don't share it with anybody.
SECRET_KEY = '^w6u&$-5r3a@r&+0fe51*lz#2)b=7&-a-qr^e_4d=932o0vsx)'

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = ('django.template.loaders.filesystem.Loader',
 'django.template.loaders.app_directories.Loader')

MIDDLEWARE_CLASSES = (
    'django.middleware.cache.UpdateCacheMiddleware',
    'django.middleware.common.CommonMiddleware',
#    'django.contrib.sessions.middleware.SessionMiddleware',
#    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.middleware.cache.FetchFromCacheMiddleware',
)

CACHE_MIDDLEWARE_SECONDS = 86400 # 30 *
CACHE_MIDDLEWARE_KEY_PREFIX = "A1"

ROOT_URLCONF = 'amview.urls'

TEMPLATE_DIRS = (
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
    os.getcwd()+"/templates/",
)

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'amview.viewer'
)


## amView configuration:

# Base path of alignments, all alignments need to be below this path
ALIGNMENT_PATH = os.getcwd()+"/examples/"

# A function to modify path names, if the public URL structure doesn't match the file system structure
MODIFY_PATH = lambda x : x

# How are the alignment and annotation files called?
# When the user specifies a directory (URL ending in "/"), these filenames
# will be checked and loaded if they're there.
ALIGNMENT_FILES = ("alignment.fa",)
ANNOTATION_FILES = ("paircoil.fa",)

## Load local settings that don't belong in revision control:
## TEMPLATE_DIRS, ALIGNMENT_PATH, ALIGNMENT_FILES, ANNOTATION_FILES
os.chdir("/srv/app/src/amview/")
execfile("local_settings.py")
