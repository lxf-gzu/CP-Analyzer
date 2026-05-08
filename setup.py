from setuptools import setup, find_packages

setup(
    name="cp-analyzer",
    version="2.0.0",
    author="Xuefei Liu",
    author_email="201307129@gznu.edu.cn",
    description="A Constant-Potential Platform for Electrocatalytic Free Energy Calculations",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/lxf-gzu/CP-Analyzer",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "cp-analyzer = cp_analyzer.main:main",   # 安装后可直接运行 cp-analyzer
        ],
    },
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "matplotlib>=3.4",
        "pandas>=1.3",
        "unidecode>=1.3",
    ],
    python_requires=">=3.8",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)