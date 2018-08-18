from setuptools import setup

setup(
    name='OmicsIntegrator',
    packages=['OmicsIntegrator'],
    package_dir={'OmicsIntegrator': 'src'},
    package_data={'OmicsIntegrator': ['annotation/final_annotation.pickle', 'src/viz.jinja', 'viz.jinja']},
    version='2.3.5',
    url='https://github.com/fraenkel-lab/OmicsIntegrator2',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    license='MIT',
    author='zfrenchee',
    author_email='alex@lenail.org',
    description='',
    install_requires=[
        "numpy",
        "pandas==0.23.4",
        "networkx==2.1",
        "pcst_fast==1.0.7",
        "python-louvain",
        "goenrich",
        "sklearn"
    ],
    entry_points={
        'console_scripts': [
            'OmicsIntegrator = src.__main__:main',
        ]
    },
)

